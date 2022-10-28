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

using namespace dealii;

/**
 * @brief Small structure used to house statistics (min,max,total,average) that are used in multiple simulations
 */

class statistics
{
public:
  statistics()
    : min(0)
    , max(0)
    , total(0)
    , average(0)
  {}

  double min;
  double max;
  double total;
  double average;
};


/**
 * @brief add_statistics_to_table_handler Add statistics to a TableHandler under the indicated variable name
 */

inline void
add_statistics_to_table_handler(const std::string variable,
                                const statistics  stats,
                                TableHandler &    table)
{
  table.add_value("Variable", variable);
  table.add_value("Min", stats.min);
  table.add_value("Max", stats.max);
  table.add_value("Total", stats.total);
  table.add_value("Average", stats.average);
};

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


template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T> &                          independent_values,
  const std::string &                             independent_column_name,
  const std::vector<std::vector<Tensor<1, dim>>> &dependent_vector,
  const std::vector<std::string> &                dependent_column_name,
  const unsigned int                              display_precision);

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
  // VOF advection equation
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


/**
 * @brief Reads a file that was built by writing a deal.II TableHandler class, and refills a TableHandler with the data in the file.
 * * @param table The table to be filled. Warning ! if the table is empty, it's content will be erased.
 *
 *    @param file The path the file that will be use to fill up the table.
 *
 *   @param delimiter The delimiter used to read the table.
 */
void
fill_table_from_file(TableHandler &    table,
                     const std::string file_name,
                     const std::string delimiter = " ");

/**
 * @brief function that read a file that was build from a dealii table and fill 2 vectors.
 * The first vector contains all the columns names and the second one contained
 * all the column data.
 * * @param map A map used to contain the data based on the columns name.
 *
 *   @param file The path the file that will be use to fill up the table.
 *
 *   @param delimiter The delimiter used to read the table.
 */
void
fill_vectors_from_file(std::map<std::string, std::vector<double>> &map,
                       const std::string                           file_name,
                       const std::string delimiter = " ");

/**
 * @brief Creates the simulation output folder
 * * @param dirname Output directory name
 */
void
create_output_folder(const std::string &dirname);



#endif
