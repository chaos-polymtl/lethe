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

#ifndef lethe_density_model_h
#define lethe_density_model_h

#include <core/physical_property_model.h>

/**
 * @brief DensityModel. Abstract class that allows to calculate the
 * density.
 */
class DensityModel : public PhysicalPropertyModel
{
public:
  /**
   * @brief Instantiates and returns a pointer to a DensityModel object by casting it to
   * the proper child class
   *
   * @param material_properties Parameters for a single fluid
   */
  static std::shared_ptr<DensityModel>
  model_cast(const Parameters::Material &material_properties);

  virtual double
  get_psi() const
  {
    return 0;
  }

  virtual double
  get_density_ref() const
  {
    return 0;
  }
};


/**
 * @brief Constant density.
 */
class DensityConstant : public DensityModel
{
public:
  /**
   * @brief Default constructor
   */
  DensityConstant(const double p_density)
    : density(p_density)
  {}

  /**
   * @brief value Calculates the density
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> & /*fields_value*/) override
  {
    return density;
  }

  /**
   * @brief vector_value Calculates the vector of density.
   * @param field_vectors Vectors of the fields on which the density may depend.
   * @param property_vector Vectors of the density values
   */
  void
  vector_value(const std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector) override
  {
    std::fill(property_vector.begin(), property_vector.end(), density);
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the density with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the density with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/,
           field /*id*/) override
  {
    return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the density with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the density with respect to the field id.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field /*id*/,
    std::vector<double> &jacobian_vector) override
  {
    std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

private:
  const double density;
};


/**
 * @brief Isothermal ideal gas density.
 */
class DensityIsothermalIdealGas : public DensityModel
{
public:
  /**
   * @brief Default constructor
   */
  DensityIsothermalIdealGas(const double p_density_ref,
                            const double p_R,
                            const double p_T)
    : density_ref(p_density_ref)
    , psi(1. / (p_R * p_T))
  {
    this->model_depends_on[field::pressure] = true;
  }

  /**
   * @brief value Calculates the density
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  double
  value(const std::map<field, double> &fields_value) override
  {
    const double pressure = fields_value.at(field::pressure);
    return density_ref + psi * pressure;
  }

  /**
   * @brief vector_value Calculates the vector of density.
   * @param field_vectors Vectors of the fields on which the density may depend.
   * @param property_vector Vectors of the density values
   */
  void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double> &property_vector) override
  {
    const std::vector<double> &pressure = field_vectors.at(field::pressure);
    for (unsigned int i = 0; i < property_vector.size(); ++i)
      property_vector[i] = density_ref + psi * pressure[i];
  }

  /**
   * @brief jacobian Calculates the jacobian (the partial derivative) of the density with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated.
   * @return value of the partial derivative of the density with respect to the field.
   */
  double
  jacobian(const std::map<field, double> & /*field_values*/, field id) override
  {
    if (id == field::pressure)
      return psi;
    else
      return 0;
  }

  /**
   * @brief vector_jacobian Calculates the derivative of the density with respect to a field.
   * @param field_vectors Vector for the values of the fields used to evaluate the property.
   * @param id Identifier of the field with respect to which a derivative should be calculated.
   * @param jacobian vector of the value of the derivative of the density with respect to the field id.
   */
  void
  vector_jacobian(
    const std::map<field, std::vector<double>> & /*field_vectors*/,
    const field          id,
    std::vector<double> &jacobian_vector) override
  {
    if (id == field::pressure)
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), psi);
    else
      std::fill(jacobian_vector.begin(), jacobian_vector.end(), 0);
  }

  /**
   * @brief get_psi Returns the value of the compressibility factor used in the density model
   * @return isothermal ideal gas compressibility factor
   */
  double
  get_psi() const override
  {
    return psi;
  }

  /**
   * @brief get_density_ref Returns the value of the reference state density used in the density model
   * @return isothermal ideal gas reference state density
   */
  double
  get_density_ref() const override
  {
    return density_ref;
  }

private:
  const double density_ref;
  const double psi;
};

#endif
