// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

using namespace dealii;

/**
 * @brief Construct oversampled quadrature
 *
 * @param[in] quadrature Quadrature for local cell operations
 *
 * @return Quadrature oversampled
 */
template <int dim>
Quadrature<dim>
construct_quadrature(const Quadrature<dim> &quad)
{
  const double oversampling_factor = 2.0; // make parameter

  for (unsigned int i = 1; i <= 10; ++i)
    if (quad == QGauss<dim>(i))
      return QGauss<dim>(i * oversampling_factor);

  AssertThrow(false, ExcNotImplemented());

  return quad;
}

/**
 * @brief Compute inner product
 *
 * @param[in] grad Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-0 tensor
 */
template <int dim, typename Number>
Number
contract(const Tensor<1, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

/**
 * @brief Compute inner product
 *
 * @param[in] grad Rank-2 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-1 tensor
 */
template <int dim, typename Number>
Tensor<1, dim, Number>
contract(const Tensor<2, dim, Number> &grad,
         const Tensor<1, dim, Number> &normal)
{
  return grad * normal;
}

/**
 * @brief Compute outer product for the multicomponent case
 *
 * @param[in] grad Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-0 tensor
 *
 */
template <int n_components, int dim, typename Number>
Tensor<1, n_components, Number>
contract(const Tensor<1, n_components, Tensor<1, dim, Number>> &grad,
         const Tensor<1, dim, Number>                          &normal)
{
  Tensor<1, n_components, Number> result;

  for (int r = 0; r < n_components; ++r)
    result[r] = grad[r] * normal;

  return result;
}

/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-0 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-1 tensor
 */
template <int dim, typename Number>
Tensor<1, dim, Number>
outer(const Number &value, const Tensor<1, dim, Number> &normal)
{
  return value * normal;
}


/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-1 tensor
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-2 tensor
 */
template <int dim, typename Number>
Tensor<2, dim, Number>
outer(const Tensor<1, dim, Number> &value, const Tensor<1, dim, Number> &normal)
{
  Tensor<2, dim, Number> result;

  for (int c = 0; c < dim; ++c)
    result[c] = value[c] * normal;

  return result;
}

/**
 * @brief Compute outer product
 *
 * @param[in] value Rank-1 tensor for n_components
 * @param[in] normal Rank-1 tensor
 *
 * @return Rank-2 tensor for n_components
 */
template <int n_components, int dim, typename Number>
Tensor<1, n_components, Tensor<1, dim, Number>>
outer(const Tensor<1, n_components, Number> &value,
      const Tensor<1, dim, Number>          &normal)
{
  Tensor<1, n_components, Tensor<1, dim, Number>> result;

  for (int c = 0; c < n_components; ++c)
    result[c] = value[c] * normal;

  return result;
}

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-1 tensor within rank-1 tensor where result is
 * stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
template <int dim, int dim_, typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<1, dim_, Tensor<1, dim, Number>> &v_gradient,
                        const Tensor<2, dim, Number>            &u_gradient,
                        const Number                            &factor)
{
  for (int d = 0; d < dim; ++d)
    v_gradient[d][d] += u_gradient[d][d] * factor;

  for (int e = 0; e < dim; ++e)
    for (int d = e + 1; d < dim; ++d)
      {
        const auto tmp = (u_gradient[d][e] + u_gradient[e][d]) * (factor * 0.5);
        v_gradient[d][e] += tmp;
        v_gradient[e][d] += tmp;
      }
}

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-2 tensor where result is stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
template <int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<2, dim, Number>       &v_gradient,
                        const Tensor<2, dim, Number> &u_gradient,
                        const Number                 &factor)
{
  for (int d = 0; d < dim; ++d)
    v_gradient[d][d] += u_gradient[d][d] * factor;

  for (int e = 0; e < dim; ++e)
    for (int d = e + 1; d < dim; ++d)
      {
        const auto tmp = (u_gradient[d][e] + u_gradient[e][d]) * (factor * 0.5);
        v_gradient[d][e] += tmp;
        v_gradient[e][d] += tmp;
      }
}

/**
 * @brief Compute scalar product
 *
 * @param[in, out] v_gradient Rank-2 tensor where result is stored
 * @param[in] u_gradient Rank-2 tensor
 * @param[in] factor Scalar factor
 */
template <typename Number>
inline DEAL_II_ALWAYS_INLINE void
symm_scalar_product_add(Tensor<1, 1, Number>       &v_gradient,
                        const Tensor<1, 1, Number> &u_gradient,
                        const Number               &factor)
{
  v_gradient[0] += u_gradient[0] * factor;
}


template <typename T>
class BufferRW
{
public:
  BufferRW(T *ptr, const unsigned int offset)
    : ptr(ptr ? (ptr + offset) : nullptr)
  {}

  void
  write(const T &in)
  {
    ptr[0] = in;
    ptr += 1;
  }

  template <int dim>
  void
  write(const Tensor<1, dim, T> &in)
  {
    for (int i = 0; i < dim; ++i)
      ptr[i] = in[i];

    ptr += dim;
  }

  template <typename T0>
  T0
  read() const
  {
    T0 result = {};

    if (ptr)
      read(result);

    return result;
  }

private:
  mutable T *ptr;

  template <int dim>
  void
  read(Tensor<1, dim, T> &out) const
  {
    for (int i = 0; i < dim; ++i)
      out[i] = ptr[i];

    ptr += dim;
  }

  void
  read(T &out) const
  {
    out = ptr[0];
    ptr += 1;
  }
};
