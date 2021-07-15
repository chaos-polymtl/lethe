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

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2020 -
 */

#ifndef lethe_utilities_h
#define lethe_utilities_h

#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping_manifold.h>
#include <deal.II/fe/mapping_q1.h>

using namespace dealii;

template <int dim>
typename DoFHandler<dim>::active_cell_iterator
find_cell_around_point_with_tree(const DoFHandler<dim> &dof_handler,
                                 Point<dim>             point);


// Utility function to create tables from vectors of scalars/tensors as
// dependent or independent variables.



/**
 * @brief Generate a table from a vector of scalar and a vector of tensor<1,dim>
 *
 * @param independent_values A vector of scalar values that serve as the independent
 * variable. For example time.
 *
 * @param independent_column_name The name of the independent variable
 *
 * @param dependent_vector. A vector of Tensor<1,dim> which are the dependent variable. For example force.
 *
 * @param dependent_column_names. A vector of string which are the label of dependent tensor
 *
 * @param display_precision. An integer which indicate the precision at which the tables are written
 *
 */
template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T> &             independent_values,
  const std::string &                independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
  const unsigned int                 display_precision);

/**
 * @brief Generate a table from a vector of tensor<1,dim> and a vector of tensor<1,dim>
 *
 * @param independent_values A vector of Tensor<1,dim> that serve as the independent
 * variable. For example position.
 *
 * @param independent_column_name A vector of string which are the label of the independent tensor.
 *
 * @param dependent_vector. A vector of Tensor<1,dim> which are the dependent variable. For example force.
 *
 * @param dependent_column_names. A vector of string which are the label of dependent tensor.
 *
 * @param display_precision. An integer which indicate the precision at which the tables are written.
 *
 */
template <int dim>
TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string> &   independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
  const unsigned int                 display_precision);


/**
 * @brief Generate a table from a vector of tensor<1,dim> and a vector of tensor<1,dim>.
 *
 * @param independent_values A vector of Tensor<1,dim> that serve as the independent
 * variable. For example position.
 *
 * @param independent_column_name The name of the independent variable.
 *
 * @param dependent_vector. A vector of scalar which are the dependent variable. For example energy.
 *
 * @param dependent_column_names. The label of the dependent scalar.
 *
 * @param display_precision. An integer which indicate the precision at which the tables are written.
 *
 */
template <int dim>
TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string> &   independent_column_name,
  const std::vector<double> &        dependent_values,
  const std::string &                dependent_column_name,
  const unsigned int                 display_precision);



/**
 * @brief Used in multiphasic simulations to calculates the equivalent properties for a given phase. Method called in quadrature points loop
 *
 * @param phase Phase value for the given quadrature point
 *
 * @param property0 Property value for the fluid with index 0 (fluid for phase = 0)
 *
 * @param property1 Property value for the fluid with index 1 (fluid for phase = 1)
 */
inline double
calculate_point_property(const double phase,
                         const double property0,
                         const double property1)
{
  double property_eq = phase * property1 + (1 - phase) * property0;

  // Limit parameters value (patch)
  // TODO see if necessary after compression term is added in the
  // free_surface advection equation
  const double property_min = std::min(property0, property1);
  const double property_max = std::max(property0, property1);
  if (property_eq < property_min)
    {
      property_eq = property_min;
    }
  if (property_eq > property_max)
    {
      property_eq = property_max;
    }

  return property_eq;
}

#endif
