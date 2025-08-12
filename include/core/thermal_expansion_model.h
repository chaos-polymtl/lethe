// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_thermal_expansion_model_h
#define lethe_thermal_expansion_model_h

#include <core/phase_change.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * thermal expansion coefficient (with unit of inverse temperature)
 */
class ThermalExpansionModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ThermalExpansionModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ThermalExpansionModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant thermal expansion coefficient.
 */
class ConstantThermalExpansion : public ThermalExpansionModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantThermalExpansion(const double p_thermal_expansion)
    : thermal_expansion(p_thermal_expansion)
  {}

  /**
   * @brief value Calculates the thermale expansion
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return thermal_expansion;
  };

  /**
   * @brief vector_value Calculates the vector of thermal expansion.
   * @param field_vectors Vectors of the fields on which the thermal expansion may depend.
   * @param property_vector Vectors of the thermal expansion values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(),
              property_vector.end(),
              thermal_expansion);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the thermal expansion with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the thermal expansion with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculates the derivative of the thermal expansion with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the thermal expansion with respect to the field id.
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
  const double thermal_expansion;
};

/**
 * @brief Implements a phase-dependent thermal expansion coefficient
 */
class ThermalExpansionPhaseChange : public ThermalExpansionModel
{
public:
  /**
   * @brief Default constructor
   */
  ThermalExpansionPhaseChange(const Parameters::PhaseChange phase_change_params)
    : p_phase_change_params(phase_change_params)
  {
    this->model_depends_on[field::temperature] = true;
  }

  /**
   * @brief value Calculates the value of the thermal expansion coefficient
   * @param fields_value Value of the various field on which the thermal expansion coefficient depends.
   * @return value of the thermal expansion coefficient calculated with the fields_value.
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    Assert(fields_value.find(field::temperature) != fields_value.end(),
           PhysicialPropertyModelFieldUndefined("ThermalExpansionPhaseChange",
                                                "temperature"));
    if (fields_value.at(field::temperature) > p_phase_change_params.T_liquidus)
      return p_phase_change_params.thermal_expansion_l;
    else
      return p_phase_change_params.thermal_expansion_s;
  };

  /**
   * @brief vector_value Calculates the vector value of thermal expansion coefficients
   * @param field_vectors Vector of properties on which the thermal expansion coefficients depend
   * @param property_vector Values of the thermal expansion coefficients
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    Assert(field_vectors.find(field::temperature) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("ThermalExpansionPhaseChange",
                                                "temperature"));
    const std::vector<double> &T = field_vectors.at(field::temperature);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      {
        // Thermal expansion of solid phase
        if (T[i] > p_phase_change_params.T_liquidus)
          property_vector[i] = p_phase_change_params.thermal_expansion_l;
        else
          property_vector[i] = p_phase_change_params.thermal_expansion_s;
        // Thermal expansion of liquid phase
      }
  };

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the thermal expansion with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the thermal expansion with respect to the field.
   */

  double
  jacobian(const std::map<field, double> &field_values, field id) override
  {
    if (id == field::temperature)
      {
        Assert(field_values.find(field::temperature) != field_values.end(),
               PhysicialPropertyModelFieldUndefined(
                 "ThermalExpansionPhaseChange", "temperature"));
        return this->numerical_jacobian(field_values, field::temperature);
      }
    else
      return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the thermal expansion with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the thermal expansion with respect to the field id
   */

  void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) override
  {
    Assert(field_vectors.find(field::temperature) != field_vectors.end(),
           PhysicialPropertyModelFieldUndefined("ThermalExpansionPhaseChange",
                                                "temperature"));
    vector_numerical_jacobian(field_vectors, id, jacobian_vector);
  };

private:
  Parameters::PhaseChange p_phase_change_params;
};

#endif
