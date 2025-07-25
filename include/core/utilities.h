// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_utilities_h
#define lethe_utilities_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/revision.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <map>
#include <regex>
#include <string>
#include <tuple>
#include <vector>

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
 * @brief Generate a table from a vector of scalars (independent variable) and
 * a vector of vectors of scalars (dependent variables).
 *
 * @tparam T Scalar type of independent variables.
 *
 * @param[in] independent_values Vector of scalar values that serve as the
 * independent variable (e.g., time).
 *
 * @param[in] independent_column_name Label of the independent variable.
 *
 * @param[in] dependent_vector Vector of vectors of scalar values containing
 * dependant variables values.
 *
 * @param[in] dependent_column_names Vector of strings representing the labels
 * of dependent variables.
 *
 * @param[in] display_precision Integer indicating the precision at which the
 * table is written.
 *
 * @param[in] display_scientific_notation Boolean indicating if the values
 * should be displayed in scientific notation (true) or not (false).
 *
 * @return Table with the independent variable values followed by the dependent
 * variables values.
 */
template <typename T>
TableHandler
make_table_scalars_vectors(
  const std::vector<T>                   &independent_values,
  const std::string                      &independent_column_name,
  const std::vector<std::vector<double>> &dependent_vector,
  const std::vector<std::string>         &dependent_column_names,
  const unsigned int                      display_precision,
  const bool                              display_scientific_notation = false);


/**
 * @brief Generate a table from a vector of scalars (independent variable) and
 * a vector of Tensor<1,dim> (dependent variables).
 *
 * @tparam T Scalar type of independent variables.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @param[in] independent_values Vector of scalar values that serve as the
 * independent variable (e.g., time).
 *
 * @param[in] independent_column_name Label of the independent variable.
 *
 * @param[in] dependent_vector Vector of Tensor<1,dim> containing dependent
 * variables values (e.g., force).
 *
 * @param[in] dependent_column_names Vector of strings representing the labels
 * of dependent variables.
 *
 * @param[in] display_precision Integer indicating the precision at which the
 * table is written.
 *
 * @param[in] display_scientific_notation Boolean indicating if the values
 * should be displayed in scientific notation (true) or not (false).
 *
 * @return Table with the independent variable values followed by the dependent
 * variables values.
 */
template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T>              &independent_values,
  const std::string                 &independent_column_name,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_names,
  const unsigned int                 display_precision,
  const bool                         display_scientific_notation = false);

/**
 * @brief Generate a table from a vector of scalar (independent variable) and a
 * vector of vectors of Tensor<1,dim> (dependent variables).
 *
 * @tparam T Scalar type of independent variables.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @param[in] independent_values Vector of scalar values that serve as the
 * independent variable (e.g., time).
 *
 * @param[in] independent_column_name Label of the independent variable.
 *
 * @param[in] dependent_vector Vector of vectors of Tensor<1,dim> containing
 * dependent variables values (e.g., force).
 *
 * @param[in] dependent_column_names Vector of strings representing the labels
 * of dependent variables.
 *
 * @param[in] display_precision Integer indicating the precision at which the
 * table is written.
 *
 * @param[in] display_scientific_notation Boolean indicating if the values
 * should be displayed in scientific notation (true) or not (false).
 *
 * @return Table with the independent variable values followed by the dependent
 * variables values.
 */
template <int dim, typename T>
TableHandler
make_table_scalars_tensors(
  const std::vector<T>                           &independent_values,
  const std::string                              &independent_column_name,
  const std::vector<std::vector<Tensor<1, dim>>> &dependent_vector,
  const std::vector<std::string>                 &dependent_column_names,
  const unsigned int                              display_precision,
  const bool display_scientific_notation = false);

/**
 * @brief Generate a table from a vector of Tensor<1,dim> (independent
 * variables) and a vector of Tensor<1,dim> (dependent variables).
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @param[in] independent_vector Vector of Tensor<1,dim> that serves as the
 * independent variable (e.g., position).
 *
 * @param[in] independent_column_names Vector of strings representing labels of
 * the independent tensor.
 *
 * @param[in] dependent_vector Vector of vectors of Tensor<1,dim> containing
 * dependent variables values (e.g., force).
 *
 * @param[in] dependent_column_names Vector of strings representing the labels
 * of dependent variables.
 *
 * @param[in] display_precision Integer indicating the precision at which the
 * table is written.
 *
 * @param[in] display_scientific_notation Boolean indicating if the values
 * should be displayed in scientific notation (true) or not (false).
 *
 * @return Table with the independent variables values followed by the dependent
 * variables values.
 */
template <int dim>
TableHandler
make_table_tensors_tensors(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string>    &independent_column_names,
  const std::vector<Tensor<1, dim>> &dependent_vector,
  const std::vector<std::string>    &dependent_column_names,
  const unsigned int                 display_precision,
  const bool                         display_scientific_notation = false);


/**
 * @brief Generate a table from a vector of Tensor<1,dim> (independent
 * variables) and a vector of scalars (dependent variable).
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 *
 * @param[in] independent_vector Vector of Tensor<1,dim> that serves as the
 * independent variable (e.g., position).
 *
 * @param[in] independent_column_names Vector of strings representing labels of
 * the independent tensor.
 *
 * @param[in] dependent_values Vector of of doubles containing
 * dependent variables values (e.g., force).
 *
 * @param[in] dependent_column_name Label of the dependent variable.
 *
 * @param[in] display_precision Integer indicating the precision at which the
 * table is written.
 *
 * @param[in] display_scientific_notation Boolean indicating if the values
 * should be displayed in scientific notation (true) or not (false).
 *
 * @return Table with the independent variables values followed by the dependent
 * variable values.
 */
template <int dim>
TableHandler
make_table_tensors_scalars(
  const std::vector<Tensor<1, dim>> &independent_vector,
  const std::vector<std::string>    &independent_column_names,
  const std::vector<double>         &dependent_values,
  const std::string                 &dependent_column_name,
  const unsigned int                 display_precision,
  const bool                         display_scientific_notation = false);


/**
 * @brief Calculate the equivalent properties for a given phase. Method called
 * in quadrature points loops in VOF simulations.
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
 * @brief Clip a property between a lower and a higher value
 * @tparam T argument's type
 * @param[in,out] lower Lower admissible value
 * @param[in,out] upper Upper admissible value
 * @param[in,out] n
 * @return Clipped variable that is not below the lower limit and not above the upper limit
 *
 */
template <typename T>
T
clip(const T &n, const T &lower, const T &upper)
{
  return std::max(lower, std::min(n, upper));
}

/**
 * @brief Calculate the equivalent properties for a given phase. Method called
 * in quadrature points loops in Cahn-Hilliard simulations.
 *
 * @param phase_cahn_hilliard Phase value for the given quadrature point
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
 * @param table The table to be filled. Warning ! if the table is empty, it's content will be erased.
 *
 * @param file_name The path the file that will be use to fill up the table.
 *
 * @param delimiter The delimiter used to read the table.
 */
void
fill_table_from_file(TableHandler      &table,
                     const std::string &file_name,
                     const std::string &delimiter = " ");

/**
 * @brief function that read a file that was build from a dealii table and fill 2 vectors.
 * The first vector contains all the columns names and the second one contained
 * all the column data.
 * @param map A map used to contain the data based on the columns name.
 *
 * @param file_name The path the file that will be use to fill up the table.
 *
 * @param delimiter The delimiter used to read the table.
 */
void
fill_vectors_from_file(std::map<std::string, std::vector<double>> &map,
                       const std::string                          &file_name,
                       const std::string &delimiter = " ");

/**
 * @brief Function that read a file that was build from a dealii table and create a map with the key being the column name and the variable the vectors of data.
 * @param map Container to be filled with the information in the dealii table.
 *
 * @param file_name The path the file that will be use to fill up the table.
 *
 * @param delimiter The delimiter used to read the table.
 */
void
fill_string_vectors_from_file(
  std::map<std::string, std::vector<std::string>> &map,
  const std::string                               &file_name,
  const std::string                               &delimiter = " ");

/**
 * @brief Creates the simulation output folder
 * @param dirname Output directory name
 */
void
create_output_folder(const std::string &dirname);

/**
 * @brief Delete the output folder and its content.
 *
 * @param[in] dirname Output directory name.
 */
void
delete_output_folder(const std::string &dirname);

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
 * @brief Extract the maximum value between the number of boundary conditions and the number of manifolds from the file.
This value is linked to the "number" string defined in the simulation parameter
file. It provides an estimate for the amount of parameters or manifolds and is
used to determine the size of the vectors that will store boundary conditions
and manifold data. This feature will need to be monitored extensively in the
future.
 * @param file_name The file name from which the number of boundary conditions
 * is read
 */
int
get_max_subsection_size(const std::string &file_name);

/**
 * @brief Return the tensor corresponding to the @p value_string. If the
 * dimension correspondence of the @p value_string is not equivalent to
 * @p spacedim (either 2 or 3), an exception will be thrown. The delimiter
 * separating the elements of the @p value_string is a comma (",").
 *
 * @remark The function can be use to construct Point<spacedim> objects.
 *
 * @tparam spacedim Number of spatial dimensions (2D or 3D).
 *
 * @param[in] value_string A string in the parameter file corresponding to a
 * given tensor.
 *
 * @return A Tensor<1,spacedim> corresponding to the @p value_string in the
 * parameter file.
 */
template <int spacedim>
inline Tensor<1, spacedim>
value_string_to_tensor(const std::string &value_string)
{
  std::vector<std::string> vector_of_string(
    Utilities::split_string_list(value_string));
  std::vector<double> vector_of_double =
    Utilities::string_to_double(vector_of_string);

  AssertThrow(vector_of_double.size() == 3 || vector_of_double.size() == 2,
              ExcMessage("Invalid string: " + value_string +
                         ". This should be a two or three dimensional vector "
                         "or point."));

  Tensor<1, spacedim> output_tensor;
  for (unsigned int i = 0; i < spacedim; ++i)
    output_tensor[i] = vector_of_double[i];

  return output_tensor;
}

/**
 * @brief Return the tensor corresponding to the @p value_string_0, but it can
 * also allow the usage of deprecated parameters that used to be 3 individual
 * entries instead of a list of values.
 * In the case of a single entry declaration, the delimiter separating the
 * elements of the @p value_string_0 is a comma (",").
 *
 * @remark The function can be use to construct Point<spacedim> objects.
 *
 * @tparam spacedim Number of spatial dimensions (2D or 3D).
 *
 * @param[in] value_string_0 A string in the parameter file corresponding to the
 * first component of the tensor or to the tensor itself.
 * @param[in] value_1 A double in the parameter file corresponding to the
 * second component of the tensor.
 * @param[in] value_2 A double in the parameter file corresponding to the
 * third component of the tensor. Only specify if @p spacedim = 3.
 *
 * @return A Tensor<1,spacedim> corresponding to the input parameters in the
 * parameter file.
 */
template <int spacedim>
inline Tensor<1, spacedim>
value_string_to_tensor(const std::string &value_string_0,
                       const double      &value_1,
                       const double      &value_2 = 0)
{
  std::vector<std::string> vector_of_string(
    Utilities::split_string_list(value_string_0));
  Tensor<1, spacedim> output_tensor;

  // The used parameter is a list of values
  if (vector_of_string.size() > 1)
    {
      std::vector<double> vector_of_double =
        Utilities::string_to_double(vector_of_string);
      AssertThrow(vector_of_double.size() == 3 || vector_of_double.size() == 2,
                  ExcMessage(
                    "Invalid string: " + value_string_0 +
                    ". This should be a two or three dimensional vector "
                    "or point."));
      for (unsigned int i = 0; i < vector_of_double.size(); ++i)
        output_tensor[i] = vector_of_double[i];
    }
  else // Depreciated individual entries
    {
      // Since the first parameter is the alias of the new parameter,
      // the value of the first parameter is obtained for its entry
      output_tensor[0] = Utilities::string_to_double(value_string_0);
      output_tensor[1] = value_1;
      if constexpr (spacedim == 3)
        output_tensor[2] = value_2;
    }
  return output_tensor;
}

/**
 * @brief Computes equivalent cell diameter by comparing the area to a disk (2D)
 * or the volume to a sphere (3D).
 *
 * @tparam dim Number of spatial dimensions (2D or 3D).
 *
 * @param[in] cell_measure Area (2D) or volume (3D) of the cell.
 *
 * @param[in] fe_degree Polynomial degree of the shape function.
 *
 * @return Cell diameter value.
 */
template <int dim>
inline double
compute_cell_diameter(const double cell_measure, const unsigned int fe_degree)
{
  double h;
  if constexpr (dim == 2)
    h = std::sqrt(4. * cell_measure / numbers::PI) / fe_degree;
  else if constexpr (dim == 3)
    h = std::cbrt(6. * cell_measure / numbers::PI) / fe_degree;
  else
    Assert(
      false,
      ExcMessage(std::string(
        "'dim' should have a value of either 2 or 3. Only 2D and 3D simulations "
        "are supported.")));
  return h;
}

/**
 * @brief Computes the area (2D) or volume (3D) of the cell by integrating 1
 * over the cell, by summing JxW values (quadrature weights) returned by the
 * FEValues object.
 *
 * @tparam dim Number of spatial dimensions (2D or 3D).
 *
 * @param[in] JxW_values Vector of mapped quadrature weights.
 *
 * @return Area (2D) volume (3D) of the cell.
 */
inline double
compute_cell_measure_with_JxW(const std::vector<double> &JxW_values)
{
  double cell_measure = 0;
  for (const double &JxW : JxW_values)
    cell_measure += JxW;
  return cell_measure;
}

/**
 * @brief Identify minimum cell size. The cell size corresponds to the diameter
 * of a disk (2D) or sphere (3D) of equivalent area (2D) or volume (3D).
 *
 * @tparam dim Number of spatial dimensions (2D or 3D).
 *
 * @param[in] mapping Interface for transformation from the reference (unit)
 * cell to the real cell.
 *
 * @param[in] dof_handler Describes the layout of DoFs and the type of FE.
 *
 * @param[in] cell_quadrature Local cell quadrature information.
 *
 * @param[in] mpi_communicator Allows communication between cores.
 *
 * @return Smallest cell size value.
 */
template <int dim>
inline double
identify_minimum_cell_size(const Mapping<dim>    &mapping,
                           const DoFHandler<dim> &dof_handler,
                           const Quadrature<dim> &cell_quadrature,
                           const MPI_Comm        &mpi_communicator)
{
  // Initialize FEValues for interface algebraic reinitialization
  FEValues<dim> fe_values(mapping,
                          dof_handler.get_fe(),
                          cell_quadrature,
                          update_JxW_values);

  // Initialize cell diameter
  double h = DBL_MAX;

  // Element degree
  const double degree = double(fe_values.get_fe().degree);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          // Compute cell diameter
          double cell_measure =
            compute_cell_measure_with_JxW(fe_values.get_JxW_values());
          double h_local = compute_cell_diameter<dim>(cell_measure, degree);

          // Update cell diameter to minimum value
          h = std::min(h, h_local);
        }
    }

  // Get the minimum between all processes
  h = Utilities::MPI::min(h, mpi_communicator);

  return h;
}

/**
 * @brief Calculate the penalty factor for the SIPG method for the case when two cells are of different sizes. A harmonic mean is used to calculate an average cell extent.
 *
 * @param[in] fe_degree Finite element degree
 * @param[in] cell_extent_here Distance of the here cell. The ratio between the
 * volume of the cell and the area of the face is generally used as a metric for
 * this distance.
 * @param[in] cell_extent_there Distance of the there cell. The ratio between
 * the volume of the cell and the area of the face is generally used as a metric
 * for this distance.
 * @return double
 */

inline double
get_penalty_factor(const unsigned int fe_degree,
                   const double       cell_extent_here,
                   const double       cell_extent_there)
{
  return (fe_degree + 1) * (fe_degree + 1.) * 0.5 *
         (1. / cell_extent_here + 1. / cell_extent_there);
}

/**
 * @brief Extract the keys from map object
 *
 * @tparam key Key of the map
 * @tparam val Value of the map. This template parameter is not used.
 * @param[in] map Map object from which the keys are extracted
 * @return keys Vector of keys from the map
 */

template <typename key, typename val>
std::vector<key>
extract_keys_from_map(const std::map<key, val> &map)
{
  std::vector<key> keys;
  for (const auto &pair : map)
    keys.push_back(pair.first);
  return keys;
}

/**
 * @brief Extract the values from map object
 *
 * @tparam key Key of the map. This template parameter is not used.
 * @tparam val Value of the map.
 * @param[in] map Map object from which the keys are extracted
 * @return values Vector of values from the map
 */

template <typename key, typename val>
std::vector<val>
extract_values_from_map(const std::map<key, val> &map)
{
  std::vector<val> values;
  for (const auto &pair : map)
    values.push_back(pair.second);
  return values;
}

/**
 * @brief Concatenate words which are passed from the command line to a main function
 * into a single string in which the words are seperated by space.
 *
 * @param[in] argc number of arguments in C style.
 * @param[in] argv arguments themselves in C style.
 */

std::string
concatenate_strings(const int argc, char **argv);

/**
 * @brief Print the information of the deal.II and Lethe branches
 * used to run an application in the "git describe" format.
 *
 * @param[in] pcout Parallel output stream
 */
void
print_version_info(const ConditionalOStream &pcout);

/**
 * @brief Parse the arguments given to the application
 *
 * @param[in] argc Number of arguments
 * @param[in] argv Array of strings representing the arguments
 * @return Tuple containing a map with the command line flags
 * and a boolean set to 1 for the specific flag; and a vector of
 * strings with the rest of the arguments that are not flags; in Lethe,
 * this corresponds to the parameter file given to the application.
 */
std::tuple<std::map<std::string, bool>, std::vector<std::string>>
parse_args(int argc, char **argv);

/**
 * @brief Print the parameters given by the parameter file of the application
 *
 * @param[in] pcout Parallel output stream
 * @param[in] prm Object containing the parameters parsed from the parameter
 * file
 */
void
print_parameters_to_output_file(const ConditionalOStream &pcout,
                                const ParameterHandler   &prm,
                                const std::string        &file_name);

/**
 * @brief Delete vtu and pvd files
 */
void
delete_vtu_and_pvd_files(const std::string &output_path);

/**
 * @brief Converts point to angle (in radians)
 * @param[in] point Point in cartesian coordinates (x-y plane)
 * @param[in] center_of_rotation Center of rotation
 *
 * @return Angle between point and the x-axis
 */
template <int dim>
inline double
point_to_angle(const Point<dim> &point,
               const Point<dim> &center_of_rotation = Point<dim>())
{
  return std::fmod(std::atan2(point[1] - center_of_rotation[1],
                              point[0] - center_of_rotation[0]) +
                     2 * numbers::PI,
                   2 * numbers::PI);
}

/**
 * @brief Converts radius to point in cartesian coordinates (in the x-y plane)
 * @param[in] radius Radial distance
 * @param[in] angle Angle (in radians)
 *
 * @return point Point cartesian coordinates
 */
template <int dim>
inline Point<dim>
radius_to_point(const double radius, const double rad)
{
  Point<dim> point;

  point[0] = radius * std::cos(rad);
  point[1] = radius * std::sin(rad);

  return point;
}



#endif
