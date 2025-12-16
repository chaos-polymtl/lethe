// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_electromagnetic_permeability_model_h
#define lethe_electromagnetic_permeability_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * electromagnetic permeability on each quadrature point.
 */
class ElectroMagneticPermeabilityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ElectromagneticsPermeabilityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ElectroMagneticPermeabilityModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant electromagnetic permeability.
 */
class ConstantElectroMagneticPermeability
  : public ElectroMagneticPermeabilityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantElectroMagneticPermeability(
    const double p_electromagnetic_permeability)
    : electromagnetic_permeability(p_electromagnetic_permeability)
  {}

  /**
   * @brief value Calculates the value the electromagnetic permeability
   * @param fields_value Value of the various field on which the electromagnetic permeability depends.
   * @return value of the electromagnetic permeability calculated with the fields_value.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return electromagnetic_permeability;
  };

  /**
   * @brief vector_value Calculates the vector value of electromagnetic conductivities
   * @param field_vectors Vector of properties on which the electromagnetic conductivities depend
   * @param property_vector Values of the electromagnetic conductivities
   */
  void
  vector_value(
    [[maybe_unused]] const std::map<field, std::vector<double>> &field_vectors,
    std::vector<double> &property_vector) override
  {
    property_vector.assign(property_vector.size(),
                           electromagnetic_permeability);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the electromagnetic permeability with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the electromagnetic permeability with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the electromagnetic permeability with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the electromagnetic permeability with respect to the field id
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
  const double electromagnetic_permeability;
};

#endif
