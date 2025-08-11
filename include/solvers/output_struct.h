//
// Created by victor on 2025-08-11.
//

#ifndef LETHE_OUTPUT_STRUCT_H
#define LETHE_OUTPUT_STRUCT_H

/**
 * @brief Struct containing information about the solution output.
 * It is used to pass all the information required upon calling
 * add_data_vector() for the DataOut instance without losing track of its
 * attributes such as name and data component interpretation.
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
    : dof_handler(nullptr)
    , solution(solution)
    , data_postprocessor(data_postprocessor)
    , is_postprocessor(true)
  {}
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
  OutputStruct(
    const DoFHandler<dim>          &dof_handler,
    const VectorType               &solution,
    const std::vector<std::string> &solution_names,
    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      &data_component_interpretation)
    : dof_handler(dof_handler)
    , solution(solution)
    , solution_names(solution_names)
    , data_component_interpretation(data_component_interpretation)
    , is_postprocessor(false)
    , data_postprocessor(nullptr)
  {
    const DoFHandler<dim> &dof_handler_ref =
      dynamic_cast<const DoFHandler<dim> &>(dof_handler);
  }
  const DoFHandler<dim>   &dof_handler;
  const VectorType        &solution;
  std::vector<std::string> solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
                                                 data_component_interpretation;
  const std::shared_ptr<DataPostprocessor<dim>> &data_postprocessor;
  bool                                           is_postprocessor;
};

#endif // LETHE_OUTPUT_STRUCT_H
