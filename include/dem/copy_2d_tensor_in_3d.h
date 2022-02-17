/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <dem/dem_solver_parameters.h>

using namespace dealii;

#ifndef copy_2d_tensor_in_3d_h
#  define copy_2d_tensor_in_3d_h

/**
 * Copies a two-dimensional tensor in a three-dimensional tensor. The third
 * element of the three-dimensional tensor (in z direction) is set to zero.
 *
 * @param tensor_2d Two-dimensional input tensor
 * @return tensor_3d Three-dimensional output tensor
 *
 */
inline Tensor<1, 3>
copy_2d_tensor_in_3d(const Tensor<1, 2> tensor_2d)
{
  Tensor<1, 3> tensor_3d;
  tensor_3d[0] = tensor_2d[0];
  tensor_3d[1] = tensor_2d[1];
  tensor_3d[2] = 0.0;

  return tensor_3d;
}
inline Tensor<1, 3>
copy_2d_tensor_in_3d(const Tensor<1, 3> tensor_3d)
{
  return tensor_3d;
}
/**
 * Copies a three-dimensional tensor in a two-dimensional point. The third
 * element of the three-dimensional point (in z direction) is neglected.
 *
 * @param point_2d Two-dimensional input point
 * @return point_3d Three-dimensional output point
 *
 */
inline Tensor<1, 2>
copy_3d_tensor_in_2d(const Tensor<1, 3> tensor_3d)
{
  Tensor<1, 2> tensor_2d;
  tensor_2d[0] = tensor_3d[0];
  tensor_2d[1] = tensor_3d[1];

  return tensor_2d;
}
inline Tensor<1, 2>
copy_3d_tensor_in_2d(const Tensor<1, 2> tensor_3d)
{
  return tensor_3d;
}
/**
 * Copies a two-dimensional point in a three-dimensional point. The third
 * element of the three-dimensional point (in z direction) is set to zero.
 *
 * @param point_2d Two-dimensional input point
 * @return point_3d Three-dimensional output point
 *
 */
inline Point<3>
copy_2d_point_in_3d(const Point<2> point_2d)
{
  Point<3> point_3d;
  point_3d[0] = point_2d[0];
  point_3d[1] = point_2d[1];
  point_3d[2] = 0.0;

  return point_3d;
}
inline Point<3>
copy_2d_point_in_3d(const Point<3> point_3d)
{
  return point_3d;
}
/**
 * Copies a three-dimensional point in a two-dimensional point. The third
 * element of the three-dimensional point (in z direction) is neglected.
 *
 * @param point_2d Two-dimensional input point
 * @return point_3d Three-dimensional output point
 *
 */
inline Point<2>
copy_3d_point_in_2d(const Point<3> point_3d)
{
  Point<2> point_2d;
  point_2d[0] = point_3d[0];
  point_2d[1] = point_3d[1];

  return point_2d;
}
inline Point<2>
copy_3d_point_in_2d(const Point<2> point_3d)
{
  return point_3d;
}

#endif
