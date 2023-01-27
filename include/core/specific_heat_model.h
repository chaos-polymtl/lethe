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


#include <core/phase_change.h>
#include <core/physical_property_model.h>

using namespace dealii;

/**
 * @brief SpecificHeatModel. Abstract class that allows to calculate the
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
   * @param physical_properties Parameters for a single fluid
   */
  static std::shared_ptr<SpecificHeatModel>
  model_cast(const Parameters::Material &fluid_properties);
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
 * @brief PhaseChangeSpecificHeat. This models takes into account the phase change of a material
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
    this->model_depends_on[field::temperature]          = true;
    this->model_depends_on[field::previous_temperature] = true;
  }

  /**
   * @brief value Calculates the value of the phase change specific heat.
   * @param fields_value Value of the various field on which the specific heat depends.
   * @return value of the specific heat.   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    const double temperature = fields_value.at(field::temperature);
    const double previous_temperature =
      fields_value.at(field::previous_temperature);
    if (temperature > previous_temperature)
      {
        const double dT = std::max(temperature - previous_temperature, 1e-6);
        return (enthalpy(temperature) - enthalpy(temperature - dT)) / dT;
      }
    else
      {
        const double dT = std::max(previous_temperature - temperature, 1e-6);
        return (enthalpy(temperature + dT) - enthalpy(temperature)) / dT;
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
    const std::vector<double> &previous_temperature_vec =
      field_vectors.at(field::previous_temperature);

    const unsigned int n_values = temperature_vec.size();

    Assert(n_values == previous_temperature_vec.size(),
           SizeOfFields(n_values, previous_temperature_vec.size()));
    Assert(n_values == property_vector.size(),
           SizeOfFields(n_values, property_vector.size()));

    for (unsigned int i = 0; i < n_values; ++i)
      {
        const double temperature          = temperature_vec[i];
        const double previous_temperature = previous_temperature_vec[i];
        if (temperature > previous_temperature)
          {
            const double dT =
              std::max(temperature - previous_temperature, 1e-6);
            property_vector[i] =
              (enthalpy(temperature) - enthalpy(temperature - dT)) / dT;
          }
        else
          {
            const double dT =
              std::max(previous_temperature - temperature, 1e-6);
            property_vector[i] =
              (enthalpy(temperature + dT) - enthalpy(temperature)) / dT;
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
