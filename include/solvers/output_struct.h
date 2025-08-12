// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_output_struct_h
#define lethe_output_struct_h

#include <core/vector.h>

#include <solvers/postprocessors.h>


/**
 * @brief Struct containing information about the solution output.
 * It is used to pass all the information required upon calling
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
   * @param[in] solution Solution field.
   *
   * @param[in] data_postprocessor Data postprocessor.
   *
   */
  OutputStructPostprocessor(
    const VectorType                              &solution,
    const std::shared_ptr<DataPostprocessor<dim>> &data_postprocessor)
    : solution(solution)
    , data_postprocessor(data_postprocessor)
  {}
  const VectorType                              &solution;
  const std::shared_ptr<DataPostprocessor<dim>> &data_postprocessor;
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
struct OutputStructDoFHandler
{
  /**
   * @brief Constructor for when data is not related to a data_postprocessor.
   *
   * @param[in] dof_handler DoFHandler containing solution.
   *
   * @param[in] solution Solution field.
   *
   * @param[in] solution_names Vector with strings containing solution name.
   *
   * @param[in] data_component_interpretation Vector with data component
   * interpretation of the solution.
   */
  OutputStructDoFHandler(
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

template <int dim, typename VectorType>
using OutputStruct = std::variant<OutputStructPostprocessor<dim, VectorType>,
                                  OutputStructDoFHandler<dim, VectorType>>;

#endif
