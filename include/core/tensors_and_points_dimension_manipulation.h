// SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_tensors_and_points_dimension_manipulation_h
#define lethe_tensors_and_points_dimension_manipulation_h

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

using namespace dealii;

/**
 * Copies a two-dimensional tensor in a three-dimensional tensor. The third
 * element of the three-dimensional tensor (in z direction) is set to zero.
 * If the tensor given is already 3d, do nothing and return the tensor.
 *
 * @param tensor the tensor that should be transformed in 3d
 * @return tensor_3d Three-dimensional output tensor
 *
 */

template <int dim>
inline Tensor<1, 3>
tensor_nd_to_3d(const Tensor<1, dim> &tensor)
{
  Tensor<1, 3> tensor_3d;
  tensor_3d[0] = tensor[0];
  tensor_3d[1] = tensor[1];
  if constexpr (dim > 2)
    tensor_3d[2] = tensor[2];

  if constexpr (dim == 2)
    tensor_3d[2] = 0;

  return tensor_3d;
}

/**
 * Copies a three-dimensional tensor in a two-dimensional tensor. The third
 * element of the three-dimensional point (in z direction) is neglected.
 * If the vector is tensor given is two-dimensional, do nothing and return the
 * tensor.
 *
 * @param tensor Three-dimensional input tensor
 * @return tensor_2d Three-dimensional output tensor
 *
 */
template <int dim>
inline Tensor<1, 2>
tensor_nd_to_2d(const Tensor<1, dim> &tensor)
{
  Tensor<1, 2> tensor_2d;
  tensor_2d[0] = tensor[0];
  tensor_2d[1] = tensor[1];
  return tensor_2d;
}

/**
 * Copies a two-dimensional point in a three-dimensional point. The third
 * element of the three-dimensional point (in z direction) is set to zero.
 * If the point given is already 3d, do nothing and return the point.
 *
 * @param point Two-dimensional input point
 * @return point_3d Three-dimensional output point
 *
 */
template <int dim>
inline Point<3>
point_nd_to_3d(const Point<dim> &point)
{
  Point<3> point_3d;
  point_3d[0] = point[0];
  point_3d[1] = point[1];
  if constexpr (dim > 2)
    point_3d[2] = point[2];
  if constexpr (dim == 2)
    point_3d[2] = 0;

  return point_3d;
}

/**
 * Copies a three-dimensional point in a two-dimensional point. The third
 * element of the three-dimensional point (in z direction) is neglected.
 * If the point given is already two dimensional, do nothing and return the
 * point.
 *
 * @param point Three-dimensional input point
 * @return point_2d Two-dimensional output point
 *
 */
template <int dim>
inline Point<2>
point_nd_to_2d(const Point<dim> &point)
{
  Point<2> point_2d;
  point_2d[0] = point[0];
  point_2d[1] = point[1];
  return point_2d;
}

#endif
