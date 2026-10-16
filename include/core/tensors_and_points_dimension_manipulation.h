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

/**
 * @brief Compute the tangential component of a vector on a surface defined by
 * its unit normal vector, which removes the normal component of the vector.
 * In 3D, it computes:
 * \f[
 * \mathbf{t} = \mathbf{n} \times (\mathbf{v} \times \mathbf{n})
 * \f]
 * In 2D, the unit tangent \f$\boldsymbol{\tau}\f$ is the unit normal rotated by
 * 90 degrees (cross_product_2d), and it computes:
 * \f[
 * \mathbf{t} = (\mathbf{v} \cdot \boldsymbol{\tau}) \boldsymbol{\tau}
 * \f]
 * Both expressions are equal to
 * \f$\mathbf{v} - (\mathbf{v} \cdot \mathbf{n}) \mathbf{n}\f$.
 *
 * @note This operation gives the tangential trace of a field, for instance the
 * traces in \f$H^{-1/2}(\mathrm{curl})\f$, where only the tangential
 * components are continuous across interfaces.
 *
 * @tparam dim Spatial dimension
 * @tparam Number Type of the vector
 *
 * @param[in] tensor Vector whose tangential component is computed.
 * @param[in] normal Unit normal vector of the surface.
 *
 * @return The tangential component of the vector.
 */
template <int dim, typename Number>
DEAL_II_ALWAYS_INLINE inline Tensor<1, dim, Number>
tangential_component(const Tensor<1, dim, Number> &tensor,
                     const Tensor<1, dim>         &normal)
{
  static_assert(dim == 2 || dim == 3,
                "The tangential component is only defined in 2D and 3D.");
  if constexpr (dim == 2)
    {
      // In 2D, the unit tangent is the unit normal rotated by 90 degrees
      const Tensor<1, dim> tangent = cross_product_2d(normal);
      return (tensor * tangent) * tangent;
    }
  else
    return cross_product_3d(normal, cross_product_3d(tensor, normal));
}

#endif
