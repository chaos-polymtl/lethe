/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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
 * -------------------------------------------------------------------*/


#ifndef lethe_shape_parsing_h
#define lethe_shape_parsing_h

#include <core/shape.h>
#include <core/utilities.h>

using namespace dealii;

namespace ShapeGenerator
{
  /**
   * Initializes the shape from its type and numerical arguments
   * @param type the type of shape
   * @param shape_arguments the numerical arguments
   * @param position the position of the shape
   * @param orientation the orientation of the shape
   */
  template <int dim>
  std::shared_ptr<Shape<dim>>
  initialize_shape(const std::string         type,
                   const std::vector<double> shape_arguments,
                   const Point<dim> &        position,
                   const Tensor<1, 3> &      orientation);

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
  initialize_shape(const std::string   type,
                   const std::string   file_name,
                   const Point<dim> &  position,
                   const Tensor<1, 3> &orientation);
} // namespace ShapeGenerator

#endif // lethe_shape_parsing_h
