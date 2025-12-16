// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_electric_permittivity_model_h
#define lethe_electric_permittivity_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * electric permittivity real and imaginary part on each quadrature point.
 */
class ElectricPermittivityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ElectricPermittivityModel object by casting it to
   * the proper child class for the real part of the permittivity.
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ElectricPermittivityModel>
  model_cast_real(const Parameters::Material &material_properties);

  /**
   * @brief Instantiates and returns a pointer to a ElectricPermittivityModel object by casting it to
   * the proper child class for the imaginary part of the permittivity.
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ElectricPermittivityModel>
  model_cast_imag(const Parameters::Material &material_properties);
};


/**
 * @brief Constant electric permittivity.
 */
class ConstantElectricPermittivity : public ElectricPermittivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantElectricPermittivity(const double p_electric_permittivity)
    : electric_permittivity(p_electric_permittivity)
  {}

  /**
   * @brief value Calculates the value the electric permittivity
   * @param fields_value Value of the various field on which the electric permittivity depends.
   * @return value of the electric permittivity calculated with the fields_value.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return electric_permittivity;
  };

  /**
   * @brief vector_value Calculates the vector value of electric permittivities
   * @param field_vectors Vector of properties on which the electric permittivities depend
   * @param property_vector Values of the electric permittivities
   */
  void
  vector_value(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double> &property_vector) override
  {
    property_vector.assign(property_vector.size(), electric_permittivity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the electric permittivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the electric permittivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the electric permittivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the electric permittivity with respect to the field id
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
  const double electric_permittivity;
};

#endif
