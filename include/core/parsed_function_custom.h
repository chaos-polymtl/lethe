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


#ifndef lethe_parsed_function_custom_h
#define lethe_parsed_function_custom_h


#include <deal.II/base/config.h>

#include <core/mu_parser_internal.h>
#include <core/tensors_and_points_dimension_manipulation.h>
#include <core/utilities.h>

#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>

/**
 * @brief Wrapper based on muParser to deal with arbitrary expressions
 * @tparam n_components Number of variables (and expressions) parsed
 */
template <int n_components>
class ParsedFunctionCustom
  : public FunctionTime<double>,
    protected ::internalFunctionParserCustom::ParserImplementation<n_components>
{
public:
  /**
   * @brief Basic constructor with only the finite difference step
   * @param[in] h Step for finite difference
   */
  ParsedFunctionCustom(const double h = 1e-8);

  /**
   * @brief Declaration of the parameters for the parameter file.
   *
   * @param[in,out] prm Parameter object for parameters parsing
   */
  void
  declare_parameters(ParameterHandler &prm);

  /**
   * @brief Initialization of the muParser object using member strings.
   *
   * This functions exists to avoid code duplication.
   */
  void
  initialize();

  /**
   * @brief Parse the parameters from the parameter file.
   *
   * @param[in,out] prm Parameter object for parameters parsing
   */
  void
  parse_parameters(ParameterHandler &prm);

  /**
   * @brief Initialize the muParser object using the input strings.
   *
   * @param[in] vnames Name of all variables
   * @param[in] expression Expressions for all source terms
   * @param[in] constants_list Constants used by the expressions
   */
  void
  initialize(const std::string vnames,
             const std::string expression,
             const std::string constants_list);

  /**
   * @brief Evaluate all components at the evaluation point
   *
   * @param[in] p Evaluation point
   * @param[out] values Output values
   */
  void
  vector_value(const Tensor<1, n_components> &p,
               Tensor<1, n_components>       &values) const;

  /**
   * @brief Evaluation the value at evaluation point for one component
   *
   * @param[in] p Evaluation point
   * @param[in] component Component of interest
   * @return Value
   */
  double
  value(const Tensor<1, n_components> &p,
        const unsigned int             component = 0) const;

  /**
   * @brief Evaluate the gradient at evaluation point for one component
   *
   * @param[in] p Evaluation point
   * @param[in] comp Component of interest
   * @return Gradient
   */
  Tensor<1, n_components>
  gradient(const Tensor<1, n_components> &p, unsigned int comp = 0) const;

  /**
   * @brief Evaluation the gradient at evaluation point for all components
   * @param[in] p Evaluation point
   * @param[out] gradients Gradients
   */
  void
  vector_gradient(const Tensor<1, n_components> &p,
                  Tensor<2, n_components>       &gradients) const;

  /**
   * @brief Set the time.
   * @param[in] newtime New time value
   */
  virtual void
  set_time(const double newtime) override;

private:
  /// Autoderivative behavior h step
  double h;
  /// Tensor for autoderivative behavior
  Tensor<2, n_components> ht;

  /// muParser component: variable names
  std::string vnames;
  /// muParser component:  expressions (same number as variables)
  std::string expression;
  /// muParser component: constants
  std::string constants_list;

  /// Used to know if member variables for muParser have been initialized
  bool initialized;
};

#endif
