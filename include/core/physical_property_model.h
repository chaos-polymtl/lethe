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

#ifndef lethe_physical_property_model_h
#define lethe_physical_property_model_h

#include <core/parameters.h>

using namespace dealii;

enum field : int
{
  shear_rate,
  temperature
};

/**
 * @brief PhysicalPropertyModel. Abstract class that defines the interface for a physical property model
 * Physical property model provide an abstract interface to calculate the value
 * of a physical property or a vector of physical property value for given field
 * value. By default, the interface does not require that all (or any) fields be
 * specified. This is why a map is used to pass the dependent variable. To allow
 * for the calculation of the appropriate jacobian matrix (when that is
 * necessary) the interface also provides a jacobian function, which must
 * provide the derivative with respect to the field specified as argument.
 */
class PhysicalPropertyModel
{
public:
  /**
   * @brief value Calculates the value of a physical property.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  virtual double
  value(std::map<field, double> fields_value) = 0;

  /**
   * @brief vector_value Calculates the values of a physical property for
   * @param field_vectors
   */
  virtual void
  vector_value(std::map<field, std::vector<double>> &field_vectors,
               std::vector<double>                  &property_vector) = 0;

  /**
   * @brief jacobian Calcualtes the jacobian (the partial derivative) of the physical
   * property with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the property with respect to the field.
   */

  virtual double
  jacobian(std::map<field, double> field_values, field id) = 0;

  /**
   * @brief vector_jacobian Calculate the derivative of the property with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluated the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the property with respect to the field id
   */

  virtual void
  vector_jacobian(std::map<field, std::vector<double>> &field_vectors,
                  field                                 id,
                  std::vector<double>                  &jacobian_vector) = 0;

protected:
  /**
   * @brief numerical_jacobian Calculate the jacobian through a finite difference approximation.
   * This approach, although not preferable, is meant as a fall-back when
   * calculated the jacobian manually is too difficult.
   * @param field_vectors
   * @param field_id
   * @param jacobian_vector
   * @return
   */
  inline double
  numerical_jacobian(std::map<field, double> &field_values, field id)
  {
    double f_x   = this->value(field_values);
    auto   x_dx  = field_values;
    double dx    = std::max(1e-6 * x_dx[id], 1e-8);
    x_dx[id]     = x_dx[id] + dx;
    double f_xdx = this->value(x_dx);
    return (f_xdx - f_x) / dx;
  }

  // Map to indicate on which variables the model depends on
  std::map<field, bool> model_depends_on;
};



#endif
