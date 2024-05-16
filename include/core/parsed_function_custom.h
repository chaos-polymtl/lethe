/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

template <int n_components>
class ParsedFunctionCustom
  : public FunctionTime<double>,
    protected ::internal::FunctionParserCustom::ParserImplementation<
      n_components>
{
public:
  ParsedFunctionCustom(const double h = 1e-8);

  void
  declare_parameters(ParameterHandler &prm);

  void
  initialize();

  void
  parse_parameters(ParameterHandler &prm);

  void
  initialize(const std::string vnames,
             const std::string expression,
             const std::string constants_list);

  void
  vector_value(const Tensor<1, n_components> &p,
               Tensor<1, n_components>       &values) const;

  double
  value(const Tensor<1, n_components> &p,
        const unsigned int             component = 0) const;

  Tensor<1, n_components>
  gradient(const Tensor<1, n_components> &p, unsigned int comp = 0) const;

  void
  vector_gradient(const Tensor<1, n_components> &p,
                  Tensor<2, n_components>       &gradients) const;

  virtual void
  set_time(const double newtime) override;

private:
  // Autoderivative behavior components
  double                               h;
  std::vector<Tensor<1, n_components>> ht;

  // muParser components
  std::string vnames;
  std::string expression;
  std::string constants_list;

  bool initialized;
};

#endif
