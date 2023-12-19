/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_specific_heat_model_h
#define lethe_specific_heat_model_h


#include <core/parameters.h>
#include <core/phase_change.h>
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief Abstract class that allows to calculate the
 * specific heat on each quadrature point using the temperature of the fluid.
 * magnitude. SpecficiHeat::get_specific_heat() is a pure virtual method,
 * since it can only be calculated knowing the model for the specific
 * heat that has been specifid
 */
class SpecificHeatModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a SpecificHeatModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<SpecificHeatModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class ConstantSpecificHeat : public SpecificHeatModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantSpecificHeat(const double p_specific_heat)
    : specific_heat(p_specific_heat)
  {}

  /**
   * @brief value Calculates the value of the specific heat.
   * @param fields_value Value of the various field on which the specific heat depends.
   * @return value of the specific heat.
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return specific_heat;
  };

  /**
   * @brief vector_value Calculates the vector of specific heat.
   * @param field_vectors Vector of fields on which the specific heat may depend.
   * @param property_vector Vector of specific_heat values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(), property_vector.end(), specific_heat);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the specific heat with respect to a field
   * @param field_values Value of the various fields on which the specific heat may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the specific heat with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the specific heat with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the specific heat with respect to the field id
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double specific_heat;
};


/**
 * @brief This model takes into account the phase change of a material
 * by considering the latent heat into the specific heat over a phase change
 * interval determined by [T_solidus,T_liquidus].
 */
class PhaseChangeSpecificHeat : public SpecificHeatModel
{
public:
  /**
   * @brief Default constructor
   */
  PhaseChangeSpecificHeat(
    const Parameters::PhaseChange p_phase_change_parameters)
    : param(p_phase_change_parameters)
  {
    this->model_depends_on[field::temperature]    = true;
    this->model_depends_on[field::temperature_p1] = true;
    this->model_depends_on[field::temperature_p2] = true;
  }

  /**
   * @brief value Calculates the value of the phase change specific heat.
   * @param fields_value Value of the various field on which the specific heat depends.
   * @return value of the specific heat.   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    const double temperature    = fields_value.at(field::temperature);
    const double temperature_p1 = fields_value.at(field::temperature_p1);

    const double temperature_p2 = fields_value.at(field::temperature_p2);


    // Gather information required from the simulation control to have the time
    // histort
    std::shared_ptr<SimulationControl> simulation_control =
      get_simulation_control();

    auto method = simulation_control->get_assembly_method();

    std::vector<double> time_steps_vector =
      simulation_control->get_time_steps_vector();

    Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);


    // If change between the temperature is insufficient, backtrack to the first
    // order implementation
    if (method != Parameters::SimulationControl::TimeSteppingMethod::bdf1 &&
        std::abs(temperature - temperature_p2) < 1e-6)
      method = Parameters::SimulationControl::TimeSteppingMethod::bdf1;

    switch (method)
      {
        case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
          if (temperature > temperature_p1)
            {
              const double dT = std::max(temperature - temperature_p1, 1e-6);
              return (enthalpy(temperature) - enthalpy(temperature - dT)) / dT;
            }
          else
            {
              const double dT = std::max(temperature_p1 - temperature, 1e-6);
              return (enthalpy(temperature + dT) - enthalpy(temperature)) / dT;
            }
          break;

        case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
          {
            const double enthalpy_current = enthalpy(temperature);
            const double enthalpy_p1      = enthalpy(temperature_p1);
            const double enthalpy_p2      = enthalpy(temperature_p2);
            const double dH               = bdf_coefs[0] * enthalpy_current +
                              bdf_coefs[1] * enthalpy_p1 +
                              bdf_coefs[2] * enthalpy_p2;
            const double dT = bdf_coefs[0] * temperature +
                              bdf_coefs[1] * temperature_p1 +
                              bdf_coefs[2] * temperature_p2;
            return dH / dT;
            break;
          }

        default:
          throw(std::runtime_error(
            "Other time integration schemes are not supported by the phase change model"));
      }
  }


  /**
   * @brief vector_value Calculates the vector of specific heat.
   * @param field_vectors Vector of fields on which the specific heat may depend.
   * @param property_vector Vector of specific_heat values.
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &temperature_vec =
      field_vectors.at(field::temperature);
    const std::vector<double> &p1_temperature_vec =
      field_vectors.at(field::temperature_p1);
    const std::vector<double> &p2_temperature_vec =
      field_vectors.at(field::temperature_p2);

    const unsigned int n_values = temperature_vec.size();

    Assert(n_values == p2_temperature_vec.size(),
           SizeOfFields(n_values, p2_temperature_vec.size()));
    Assert(n_values == p1_temperature_vec.size(),
           SizeOfFields(n_values, p1_temperature_vec.size()));
    Assert(n_values == property_vector.size(),
           SizeOfFields(n_values, property_vector.size()));

    // Gather information required from the simulation control to have the time
    // history
    std::shared_ptr<SimulationControl> simulation_control =
      get_simulation_control();
    std::vector<double> time_steps_vector =
      simulation_control->get_time_steps_vector();

    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double temperature    = temperature_vec[i];
        const double temperature_p1 = p1_temperature_vec[i];
        const double temperature_p2 = p2_temperature_vec[i];

        // If change between the temperature is insufficient, backtrack to the
        // first order implementation
        auto method = simulation_control->get_assembly_method();

        if (method != Parameters::SimulationControl::TimeSteppingMethod::bdf1 &&
            std::abs(temperature - temperature_p2) < 1e-6)
          method = Parameters::SimulationControl::TimeSteppingMethod::bdf1;

        Vector<double> bdf_coefs = bdf_coefficients(method, time_steps_vector);

        switch (method)
          {
            case Parameters::SimulationControl::TimeSteppingMethod::bdf1:
              if (temperature > temperature_p1)
                {
                  const double dT =
                    std::max(temperature - temperature_p1, 1e-6);
                  property_vector[i] =
                    (enthalpy(temperature) - enthalpy(temperature - dT)) / dT;
                }
              else
                {
                  const double dT =
                    std::max(temperature_p1 - temperature, 1e-6);
                  property_vector[i] =
                    (enthalpy(temperature + dT) - enthalpy(temperature)) / dT;
                }
              break;

            case Parameters::SimulationControl::TimeSteppingMethod::bdf2:
              {
                const double enthalpy_current = enthalpy(temperature);
                const double enthalpy_p1      = enthalpy(temperature_p1);
                const double enthalpy_p2      = enthalpy(temperature_p2);
                const double dH = bdf_coefs[0] * enthalpy_current +
                                  bdf_coefs[1] * enthalpy_p1 +
                                  bdf_coefs[2] * enthalpy_p2;
                const double dT = bdf_coefs[0] * temperature +
                                  bdf_coefs[1] * temperature_p1 +
                                  bdf_coefs[2] * temperature_p2;
                property_vector[i] = dH / dT;
                break;
              }

            default:
              throw(std::runtime_error(
                "Other time integration schemes are not supported by the phase change model"));
          }
      }
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the specific heat with respect to a field
   * @param field_values Value of the various fields on which the specific heat may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the specific heat with respect to the field.
   */
  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::temperature)
      return numerical_jacobian(field_values, field::temperature);
    else
      return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the specific heat with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the specific heat with respect to the field id
   */
  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    vector_numerical_jacobian(field_vectors, id, jacobian_vector);
  };

  /**
   * @brief enthalpy Calculates the enthalpy of a phase change material for a temperature T
   * The enthalpy is defined as :
   *
   * Pure liquid
   * ------------
   * if (T>T_liquidus) : H = cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T_liquidus-T_solidus) + latent_enthalpy + cp_liquid * (T-T_liquidus)
   *
   * Liquid-solid mixture
   * -----------------
   * else if (T>T_solidus) : cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T-T_solidus) + liquid_fraction * latent_enthalpy
   *
   * Pure solid
   * ------------
   * else : cp_solid * T
   *
   * @param T temperature at which to calculate the enthalpy
   * @return Value of the enthalpy
   */

  inline double
  enthalpy(const double T)
  {
    if (T > param.T_liquidus)
      return (param.cp_s * param.T_solidus +
              0.5 * (param.cp_l + param.cp_s) *
                (param.T_liquidus - param.T_solidus) +
              param.latent_enthalpy + (T - param.T_liquidus) * param.cp_l);
    else if (T > param.T_solidus)
      {
        liquid_fraction = calculate_liquid_fraction(T, param);
        return (param.cp_s * param.T_solidus +
                0.5 * (param.cp_l + param.cp_s) * (T - param.T_solidus) +
                param.latent_enthalpy * liquid_fraction);
      }
    else
      return param.cp_s * T;
  }

private:
  const Parameters::PhaseChange param;
  double                        liquid_fraction;
};

#endif
