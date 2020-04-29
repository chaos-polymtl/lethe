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

#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#ifndef lethe_utilities_h
#  define lethe_utilities_h

using namespace dealii;

template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T> &             independent_values,
  const std::string &                independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
  const unsigned int                 display_precision);

template <int dim>
TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string> &   independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string> &   dependent_column_name,
  const unsigned int                 display_precision);

template <int dim>
TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string> &   independent_column_name,
  const std::vector<double> &        dependent_values,
  const std::string &                dependent_column_name,
  const unsigned int                 display_precision);



#endif
