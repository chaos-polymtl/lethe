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

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <regex>


using namespace dealii;

/**
 * @brief Small class that is used to store statistics (min,max,total,average) of variables that are used in simulations.
 * This small class allows us to agglomerate the statistics instead of returning
 * tuples. For example, this class is used to store the kinetic energy of the
 * particles in the DEM model.
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
                                TableHandler     &table)
{
  table.add_value("Variable", variable);
  table.add_value("Min", stats.min);
  table.add_value("Max", stats.max);
  table.add_value("Total", stats.total);
  table.add_value("Average", stats.average);
}

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
  const std::vector<T>              &independent_values,
  const std::string                 &independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_name,
  const unsigned int                 display_precision);


template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T>                           &independent_values,
  const std::string                              &independent_column_name,
  const std::vector<std::vector<Tensor<1, dim>>> &dependent_vector,
  const std::vector<std::string>                 &dependent_column_name,
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
  const std::vector<std::string>    &independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_name,
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
  const std::vector<std::string>    &independent_column_name,
  const std::vector<double>         &dependent_values,
  const std::string                 &dependent_column_name,
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
 * @brief Used in the calculate_properties_ch to retrieve the sign of the phase parameter
 * @tparam T val argument's type
 * @param val value of the variable for which we want to evaluate the sign
 * @return an integer -1 or 1 depending of the sign of the phase parameter
 */
template <typename T>
int
sgn(T val)
{
  return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}

/**
 * @brief Used in multiphasic simulations to calculates the equivalent properties for a given phase given by the Cahn-Hilliard equations. Method called in quadrature points loop
 *
 * @param phase Phase value for the given quadrature point
 *
 * @param property0 Property value for the fluid with index 0 (fluid for phase = -1)
 *
 * @param property1 Property value for the fluid with index 1 (fluid for phase = 1)
 */
inline double
calculate_point_property_cahn_hilliard(const double phase_cahn_hilliard,
                                       const double property0,
                                       const double property1)
{
  double phase = phase_cahn_hilliard;

  if (std::abs(phase_cahn_hilliard) > 1)
    {
      phase = sgn(phase_cahn_hilliard);
    }

  double property_avg  = (property0 + property1) * 0.5;
  double property_diff = (property0 - property1) * 0.5;

  double property_eq = phase * property_diff + property_avg;

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
fill_table_from_file(TableHandler     &table,
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
 * @brief function that read a file that was build from a dealii table and create a map of the with the key being the column name and the variable the vectors of data.
 * * @param map A map used to contain the data based on the columns name.
 *
 *   @param file The path the file that will be use to fill up the table.
 *
 *   @param delimiter The delimiter used to read the table.
 */
void
fill_string_vectors_from_file(
  std::map<std::string, std::vector<std::string>> &map,
  const std::string                                file_name,
  const std::string                                delimiter = " ");

/**
 * @brief Creates the simulation output folder
 * * @param dirname Output directory name
 */
void
create_output_folder(const std::string &dirname);

/**
 * @brief Prints a string and then adds a line above and below made with dashes containing as many dashes as the string has characters+1
 *
 * For example, if the string to be printed is "Tracer" the result will be:
 * -------
 * Tracer
 * -------
 *
 * @param pcout the parallel cout used to print the information
 * @param expression string that will be printed
 * @param delimiter the character used to delimit the printing. Default value is "-"
 */
inline void
announce_string(const ConditionalOStream &pcout,
                const std::string         expression,
                const char                delimiter = '-')
{
  pcout << std::string(expression.size() + 1, delimiter) << std::endl;
  pcout << expression << std::endl;
  pcout << std::string(expression.size() + 1, delimiter) << std::endl;
}

/**
 * @brief Serializes a table using boost serialization feature
 * the filename should contain the desired extension
 *
 * @param table The table to be serialized
 * @param filename The file name (including the extension) to be used
 */
inline void
serialize_table(const TableHandler &table, const std::string filename)
{
  std::ofstream                 ofile(filename);
  boost::archive::text_oarchive oa(ofile, boost::archive::no_header);
  oa << table;
}

/**
 * @brief Loads a table using boost serialization feature
 * the filename should contain the desired extension
 *
 * @param table The table to be deserialized
 * @param filename The file name (including the extension) to be used
 */
inline void
deserialize_table(TableHandler &table, const std::string filename)
{
  std::ifstream                 ifile(filename);
  boost::archive::text_iarchive ia(ifile, boost::archive::no_header);
  ia >> table;
}

/**
 * @brief  get the value of a particular parameter from the contents of the input
 * file. Return an empty string if not found. This function is used to read an
 * individual parameter for an input file. This function is adapted from ASPECT
 * and is mainly used in parsing the dim of the problem before creating the
 * whole parameter parser
 *
 * @param file_name The file name from which to read a value
 * @param parameter_name The name of the parameter
 */

std::string
get_last_value_of_parameter(const std::string &file_name,
                            const std::string &parameter_name);

/**
 * @brief Extract the dimension in which to run Lethe from the
 * the contents of the parameter file. This is something that
 * we need to do before processing the parameter file since we
 * need to know whether to use the dim=2 or dim=3 instantiation
 * of the main classes.
 *
 * @param file_name The file name from which dimension is read
 */
unsigned int
get_dimension(const std::string &file_name);

/**
 * @brief Extract the maximum number of boundary conditions from the file.
 * The number of boundary conditions is linked to the string "number" so
 * this feature will need to be monitored extensively in the future.
 * @param file_name The file name from which the number of boundary conditions
 * is read
 */
int
get_max_number_of_boundary_conditions(const std::string &file_name);

#endif
