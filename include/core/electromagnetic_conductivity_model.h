// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_em_conductivity_model_h
#define lethe_em_conductivity_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief Abstract class that allows to calculate the
 * electromagnetic conductivity on each quadrature point.
 */
class ElectroMagneticConductivityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a ElectromagneticsConductivityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a material
   */
  static std::shared_ptr<ElectroMagneticConductivityModel>
  model_cast(const Parameters::Material &material_properties);
};


/**
 * @brief Constant electromagnetic conductivity.
 */
class ConstantElectroMagneticConductivity
  : public ElectroMagneticConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ConstantElectroMagneticConductivity(const double p_em_conductivity)
    : em_conductivity(p_em_conductivity)
  {}

  /**
   * @brief value Calculates the value the electromagnetic conductivity
   * @param fields_value Value of the various field on which the electromagnetic conductivity depends.
   * @return value of the electromagnetic conductivity calculated with the fields_value.
   */
  double
  value([[maybe_unused]] const std::map<field, double> &fields_value) override
  {
    return em_conductivity;
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
    property_vector.assign(property_vector.size(), em_conductivity);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the electromagnetic conductivity with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the electromagnetic conductivity with respect to the field.
   */

  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  };

  /**
   * @brief vector_jacobian Calculate the derivative of the electromagnetic conductivity with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the electromagnetic conductivity with respect to the field id
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
  const double em_conductivity;
};

#endif
