// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_shape_parsing_h
#define lethe_shape_parsing_h

#include <core/shape.h>
#include <core/utilities.h>

using namespace dealii;

namespace ShapeGenerator
{
  /**
   * Initializes the shape from its type and arguments
   * @param type the type of shape
   * @param shape_arguments_str the raw arguments
   * @param position the position of the shape
   * @param orientation the orientation of the shape
   */
  template <int dim>
  std::shared_ptr<Shape<dim>>
  initialize_shape(const std::string  &type,
                   const std::string  &shape_arguments_str,
                   const Point<dim>   &position,
                   const Tensor<1, 3> &orientation);

  /**
   * Initializes the shape from its type and numerical arguments
   * @param type the type of shape
   * @param shape_arguments the numerical arguments
   * @param position the position of the shape
   * @param orientation the orientation of the shape
   */
  template <int dim>
  std::shared_ptr<Shape<dim>>
  initialize_shape_from_vector(const std::string         &type,
                               const std::vector<double> &shape_arguments,
                               const Point<dim>          &position,
                               const Tensor<1, 3>        &orientation);

  /**
   * Initializes the shape from its type and a text file that contains the real
   * initialization data
   * @param type the type of shape
   * @param file_name the name of the file that contains the initialization data
   * @param position the position of the shape
   * @param orientation the orientation of the shape
   */
  template <int dim>
  std::shared_ptr<Shape<dim>>
  initialize_shape_from_file(const std::string  &type,
                             const std::string  &file_name,
                             const Point<dim>   &position,
                             const Tensor<1, 3> &orientation);
} // namespace ShapeGenerator

#endif // lethe_shape_parsing_h
