// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_output_struct_h
#define lethe_output_struct_h

#include <deal.II/base/table_handler.h>

#include <deal.II/numerics/data_postprocessor.h>

using namespace dealii;

/**
 * @brief Struct containing information about the solution output.
 * It is used to pass the information required upon calling
 * add_data_vector() for the DataOut instance without losing track of its
 * attributes such as name and data component interpretation. This version of
 * the struct uses Postprocessors.
 */
template <int dim, typename VectorType>
struct OutputStructPostprocessor
{
  /**
   * @brief Constructor for when data is in data_postprocessor.
   *
   * @param[in] dof_handler DoFHandler containing solution.
   * @param[in] solution Solution field.
   * @param[in] data_postprocessor Data postprocessor.
   *
   */
  OutputStructPostprocessor(
    const DoFHandler<dim>                        &dof_handler,
    const VectorType                             &solution,
    const std::shared_ptr<DataPostprocessor<dim>> data_postprocessor)
    : dof_handler(dof_handler)
    , solution(solution)
    , data_postprocessor(data_postprocessor)
  {}
  const DoFHandler<dim>                        &dof_handler;
  const VectorType                             &solution;
  const std::shared_ptr<DataPostprocessor<dim>> data_postprocessor;
};

/**
 * @brief Struct containing information about the solution output.
 * It is used to pass all the information required upon calling
 * add_data_vector() for the DataOut instance without losing track of its
 * attributes such as name and data component interpretation. This version of
 * the struct uses the DoF handler, the present solution field, and the data
 * component interpretation vector.
 */
template <int dim, typename VectorType>
struct OutputStructSolution
{
  /**
   * @brief Constructor for when data is not related to a data_postprocessor.
   *
   * @param[in] dof_handler DoFHandler containing solution.
   * @param[in] solution Solution field.
   * @param[in] solution_names Vector with strings containing solution name.
   * @param[in] data_component_interpretation Vector with data component
   * interpretation of the solution.
   */
  OutputStructSolution(
    const DoFHandler<dim>          &dof_handler,
    const VectorType               &solution,
    const std::vector<std::string> &solution_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation)
    : dof_handler(dof_handler)
    , solution(solution)
    , solution_names(solution_names)
    , data_component_interpretation(data_component_interpretation)
  {}
  const DoFHandler<dim>         &dof_handler;
  const VectorType              &solution;
  const std::vector<std::string> solution_names;
  const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation;
};

/**
 * @brief Struct containing information about the solution output.
 * It is used to pass all the information required upon calling
 * add_data_vector() for the DataOut instance without losing track of its
 * attributes such as name and data component interpretation. This version of
 * the struct uses a cell-based field parsed as a Vector and the solution name.
 */
struct OutputStructCellVector
{
  /**
   * @brief Constructor for when data is not related to a data_postprocessor.
   *
   * @param[in] solution Solution Vector<float>, which is stored as a copy.
   * @param[in] solution_name String containing solution name.
   */
  OutputStructCellVector(const Vector<float> &solution,
                         const std::string   &solution_name)
    : solution(solution)
    , solution_name(solution_name)
  {}
  const Vector<float> solution;
  const std::string   solution_name;
};

/**
 * @brief Variant handler of the three output structs
 * (OutputStructPostprocessor, OutputStructSolution and OutputStructCellVector).
 * This is used to allow the output of both postprocessors with their respective
 * DoF handlers, and the straightforward cell data struct in a single vector of
 * OutputStruct.
 */
template <int dim, typename VectorType>
using OutputStruct = std::variant<OutputStructPostprocessor<dim, VectorType>,
                                  OutputStructSolution<dim, VectorType>,
                                  OutputStructCellVector>;

/**
 * @brief Struct containing information about the table handlers to be
 * written/read. It is used to pass a vector of references to TableHandler
 * objects that should be serialized/deserialized at the writing/reading of a
 * checkpoint by the solver used, and a vector of references to strings
 * containing the names of the corresponding files.
 */
struct OutputStructTableHandler
{
  /**
   * @brief Constructor
   *
   * @param[in] table A reference to a TableHandler object that needs to be
   * serialized/ deserialized.
   * @param[in] table_filename A string with the corresponding file name.
   */
  OutputStructTableHandler(TableHandler      &table,
                           const std::string &table_filename)
    : table(table)
    , table_filename(table_filename)
  {}
  /// Reference to the TableHandler object to be serialized/deserialized
  TableHandler &table;
  /// Name of the file associated with the TableHandler for input/output
  const std::string table_filename;
};
#endif
