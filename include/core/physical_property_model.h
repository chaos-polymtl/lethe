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
#include <core/simulation_control.h>

using namespace dealii;

DeclExceptionMsg(SimulationControlUsageByProperty,
                "The usage of the simulation control object is currently only supported by the SpecificHeat models" );

DeclException2(SizeOfFields,
               unsigned int,
               unsigned int,
               << "The number of values for a field : " << arg1
               << " is not equal to the number of values for another field "
               << arg2);

enum field : int
{
  shear_rate,
  temperature,
  previous_temperature,
  pressure,
  phase_order_cahn_hilliard
};

inline void
set_field_vector(const field                          &id,
                 const std::vector<double>            &data,
                 std::map<field, std::vector<double>> &fields)
{
  std::vector<double> &target = fields.at(id);
  size_t               sz     = target.size();
  for (size_t i = 0; i < sz; ++i)
    {
      target[i] = data[i];
    }
}

/**
 * @brief Abstract class that defines the interface for a physical property model
 * Physical property model provides an abstract interface to calculate the
 * value of a physical property or a vector of physical property value for
 * given field value. By default, the interface does not require that all (or
 * any) fields be specified. This is why a map is used to pass the dependent
 * variable. To allow for the calculation of the appropriate jacobian matrix
 * (when that is necessary) the interface also provides a jacobian function,
 * which must provide the derivative with respect to the field specified as an
 * argument.
 */
class PhysicalPropertyModel
{
public:
  /**
   * @brief PhysicalPropertyModel Default constructor. Set the model_depends_on to false for all variables.
   */
  PhysicalPropertyModel()
  {
    model_depends_on[shear_rate]                = false;
    model_depends_on[temperature]               = false;
    model_depends_on[previous_temperature]      = false;
    model_depends_on[pressure]                  = false;
    model_depends_on[phase_order_cahn_hilliard] = false;
  }

  /**
   * @brief Returns true if the PhysicalPropertyModel depends on a field, false if not.
   */

  inline bool
  depends_on(const field &id)
  {
    return model_depends_on[id];
  }


  /**
   * @brief value Calculates the value of a physical property.
   * @param fields_value Value of the various field on which the property may depend.
   * @return value of the physical property calculated with the fields_value
   */
  virtual double
  value(const std::map<field, double> &fields_value) = 0;

  /**
   * @brief vector_value Calculates the values of a physical property for multiple points
   * @param field_vectors
   */
  virtual void
  vector_value(const std::map<field, std::vector<double>> &field_vectors,
               std::vector<double>                        &property_vector) = 0;

  /**
   * @brief jacobian Calcualtes the jacobian (the partial derivative) of the physical
   * property with respect to a field
   * @param field_values Value of the various fields on which the property may depend.
   * @param id Indicator of the field with respect to which the jacobian
   * should be calculated
   * @return value of the partial derivative of the property with respect to the field.
   */

  virtual double
  jacobian(const std::map<field, double> &field_values, const field id) = 0;

  /**
   * @brief vector_jacobian Calculate the derivative of the property with respect to a field
   * @param field_vectors Vector for the values of the fields used to evaluate the property
   * @param id Identifier of the field with respect to which a derivative should be calculated
   * @param jacobian Vector of the value of the derivative of the property with respect to the field id
   */

  virtual void
  vector_jacobian(const std::map<field, std::vector<double>> &field_vectors,
                  const field                                 id,
                  std::vector<double> &jacobian_vector) = 0;

  /**
   * @brief Provides the physical property with the simulation control object ensuring
   * that it can calculate physical properties that depend on time or time-history
   *
   * @param p_simulation_control shared pointed to a SimulationControl object. A copy of this shared pointer is stored in the physical property.
   */

  void
  provide_simulation_control(std::shared_ptr<SimulationControl> &p_simulation_control)
  {
    simulation_control=p_simulation_control;
  }


protected:
  /**
   * @brief numerical_jacobian Calculates the jacobian through a forward finite difference (Euler) approximation.
   * This approach, although not preferable, is meant as a fall-back when
   * calculating the jacobian manually is too difficult.
   * @param fields_values The values of the various fields
   * @param id Field id of the field with respect to which the jacobian should be calculated
   * @return
   */
  inline double
  numerical_jacobian(const std::map<field, double> &field_values,
                     const field                    id)
  {
    double f_x   = this->value(field_values);
    auto   x_dx  = field_values;
    double dx    = std::max(1e-6 * x_dx[id], 1e-8);
    x_dx[id]     = x_dx[id] + dx;
    double f_xdx = this->value(x_dx);
    return (f_xdx - f_x) / dx;
  }

  /**
   * @brief vector_numerical_jacobian Calculates the vector of jacobian through forward finite difference (Euler) approximation.
   * This approach, although not preferable, is meant as a fall-back when
   * calculating the jacobian manually is too difficult.
   * @param field_vectors
   * @param field_id
   * @param jacobian_vector
   * @return
   */
  inline void
  vector_numerical_jacobian(
    const std::map<field, std::vector<double>> &field_vectors,
    const field                                 id,
    std::vector<double>                        &jacobian_vector)
  {
    const unsigned int n_pts = jacobian_vector.size();

    // Evaluate the properties using the values of the field vector
    std::vector<double> f_x(n_pts);
    vector_value(field_vectors, f_x);

    // Make a copy of the field vector for the field we wil perturbate
    std::map<field, std::vector<double>> perturbed_field_vectors =
      field_vectors;
    std::vector<double> &x = perturbed_field_vectors.at(id);
    std::vector<double>  dx(n_pts);

    // Perturb the field by dx
    for (unsigned int i = 0; i < n_pts; ++i)
      {
        dx[i] = std::max(1e-6 * x[i], 1e-8);
        x[i] += dx[i];
      }

    // Evaluate the properties using the perturbed value of the field vector
    std::vector<double> f_xdx(n_pts);
    vector_value(perturbed_field_vectors, f_xdx);

    // Fill jacobian
    for (unsigned int i = 0; i < n_pts; ++i)
      jacobian_vector[i] = (f_xdx[i] - f_x[i]) / dx[i];
  }

  std::shared_ptr<SimulationControl> &
  get_simulation_control()
  {
    AssertThrow(simulation_control,
    SimulationControlUsageByProperty());
    return simulation_control;
  }

  // Map to indicate on which variables the model depends on
  std::map<field, bool> model_depends_on;

private:
  // SimulationControl object. This can be used to set time-dependent
  // physical properties or properties which depend on time derivatives
  // The SimulationControl must be provided to the physical property
  // through the solver itself.
  std::shared_ptr<SimulationControl> simulation_control;
};



#endif
