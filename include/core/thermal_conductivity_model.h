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

#ifndef lethe_thermal_conductivity_model_h
#define lethe_thermal_conductivity_model_h

#include <core/parameters.h>
#include <core/physical_property_model.h>

/**
 * @brief ThermalConductivityModel. Abstract class that allows to calculate the
 * thermal conductivity on each quadrature point using the temperature of the
 * fluid.
 */
class ThermalConductivityModel : public PhysicalPropertyModel
{
public:
};


/**
 * @brief Constant thermal conductivity.
 */
class ThermalConductivityConstant : public ThermalConductivityModel
{
public:
  /**
   * @brief Default constructor
   */
  ThermalConductivityConstant(const double p_thermal_conductivity)
    : thermal_conductivity(p_thermal_conductivity)
  {}

  /**
   * @brief value Calculates the value of a physical property.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  virtual double
  value(std::map<field, double> /*fields_value*/)
  {
    return thermal_conductivity;
  };

  /**
   * @brief vector_value Calculates the values of a physical property for
   * @param field_vectors
   */
  virtual void
  vector_value(std::map<field, std::vector<double>> & /*field_vectors*/,
               std::vector<double> &property_vector)
  {
    property_vector.assign(property_vector.size(), thermal_conductivity);
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
  jacobian(std::map<field, double> /*field_values*/, field /*id*/)
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
  vector_jacobian(std::map<field, std::vector<double>> & /*field_vectors*/,
                  field /*id*/,
                  std::vector<double> &jacobian_vector)
  {
    jacobian_vector.assign(jacobian_vector.size(), 0);
  };

private:
  const double thermal_conductivity;
};

#endif
