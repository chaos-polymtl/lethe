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
};


/**
 * @brief Constant specific heat. Returns a constant specific
 * heat for a fluid
 */
class SpecificHeatConstant : public SpecificHeatModel
{
public:
  /**
   * @brief Default constructor
   */
  SpecificHeatConstant(const double p_specific_heat)
    : specific_heat(p_specific_heat)
  {}

  /**
   * @brief value Calculates the value of a physical property.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  virtual double
  value(const std::map<field, double> /*fields_value*/)
  {
    return specific_heat;
  };

  /**
   * @brief vector_value Calculates the values of a physical property for
   * @param field_vectors
   */
  virtual void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector)
  {
    property_vector.assign(property_vector.size(), specific_heat);
  }

  /**
   * @brief jacobian Calcualtes the jacobian (the partial derivative) of the physical
   * property with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the property with respect to the field.
   */

  virtual double
  jacobian(const std::map<field, double> /*field_values*/, field /*id*/)
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the property with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluated the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the property with respect to the field id
   */

  virtual void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector)
  {
    jacobian_vector.assign(jacobian_vector.size(), 0);
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
  {}

  /**
   * @brief Returns the phase change specific heat
   */
  virtual double
  value(const std::map<field, double> fields_value) override
  {
    const double temperature = fields_value.at(field::temperature);
    const double previous_temperature =
      fields_value.at(field::previous_temperature);
    if (temperature > previous_temperature)
      {
        const double dT = std::max(temperature - previous_temperature, 1e-6);
        return (enthalpy(temperature) - enthalpy(previous_temperature)) / dT;
      }
    else
      {
        const double dT = std::max(previous_temperature - temperature, 1e-6);
        return (enthalpy(previous_temperature) - enthalpy(temperature)) / dT;
      }
  }


  /**
   * @brief vector_value Calculates the values of a physical property for
   * @param field_vectors
   */
  virtual void
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
              (enthalpy(temperature) - enthalpy(previous_temperature)) / dT;
          }
        else
          {
            const double dT =
              std::max(previous_temperature - temperature, 1e-6);
            property_vector[i] =
              (enthalpy(previous_temperature) - enthalpy(temperature)) / dT;
          }
      }
  }

  /**
   * @brief jacobian Calcualtes the jacobian (the partial derivative) of the physical
   * property with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the property with respect to the field.
   */

  virtual double
  jacobian(const std::map<field, double> field_values, field id) override
  {
    if (id == field::temperature)
      return numerical_jacobian(field_values, field::temperature);
    else
      return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the property with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluated the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the property with respect to the field id
   */

  virtual void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    vector_numerical_jacobian(field_vectors, id, jacobian_vector);
  };

  /**
   * @brief solid_fraction Calculates the solid fraction of a phase change material at
   * a temperature T
   *
   * @param T temperature at which to calculate the solid fraction
   *
   */
  inline double
  liquid_fraction(const double T)
  {
    return std::min(std::max((T - param.T_solidus) /
                               (param.T_liquidus - param.T_solidus),
                             0.),
                    1.);
  }


  /**
   * @brief enthalpy Calculates the enthalpy of a phase change material for a temperature T
   * The enthalpy is defined as :
   *
   * !! Pure liquid !!
   * if (T>T_liquidus) : H = cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T_liquidus-T_solidus) + latent_enthalpy + cp_liquid * (T-T_liquidus)
   *
   * !! Liquid-solid mix !!
   * else if (T>T_solidus) : cp_solid * T_solidus + 0.5*(cp_solid+cp_liquid) *
   * (T-T_solidus) + liquid_fraction * latent_enthalpy
   *
   * !! Pure solid !!
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
        const double l_frac = liquid_fraction(T);
        return (param.cp_s * param.T_solidus +
                0.5 * (param.cp_l + param.cp_s) * (T - param.T_solidus) +
                param.latent_enthalpy * l_frac);
      }
    else
      return param.cp_s * T;
  }

private:
  const Parameters::PhaseChange param;
};


#endif
