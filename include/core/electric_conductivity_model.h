// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_electric_conductivity_model_h
#define lethe_electric_conductivity_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * electric conductivity on each quadrature point.
 */
class ElectricConductivityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ElectricConductivityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ElectricConductivityModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant electric conductivity.
 */
class ConstantElectricConductivity
  : public ElectricConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantElectricConductivity(const double p_electric_conductivity)
    : electric_conductivity(p_electric_conductivity)
  {}

  /**
   * @brief value Calculates the value the electric conductivity
   * @param fields_value Value of the various field on which the electric conductivity depends.
   * @return value of the electric conductivity calculated with the fields_value.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return electric_conductivity;
  };

  /**
   * @brief vector_value Calculates the vector value of electric conductivities
   * @param field_vectors Vector of properties on which the electric conductivities depend
   * @param property_vector Values of the electric conductivities
   */
  void
  vector_value(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double> &property_vector) override
  {
    property_vector.assign(property_vector.size(),
                           electric_conductivity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the electric conductivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the electric conductivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the electric conductivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the electric conductivity with respect to the field id
   */

  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    [[maybe_unused]] const field id,
    std::vector<double>         &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  };

private:
  const double electric_conductivity;
};

#endif
