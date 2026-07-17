// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_stabilization_h
#define lethe_stabilization_h

#include <core/vector.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <cmath>

using namespace dealii;

/**
 * @brief Calculate the stabilization parameter for the Navier-Stokes equations
 * in steady-state
 * @return Value of the stabilization parameter - tau
 *
 * @param u_mag Magnitude of the velocity
 * @param kinematic_viscosity Kinematic viscosity
 * @param h Cell size; it should be calculated using the diameter of a sphere of
 * equal volume to that of the cell.
 */
inline double
calculate_navier_stokes_gls_tau_steady(const double u_mag,
                                       const double kinematic_viscosity,
                                       const double h)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(2. * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * kinematic_viscosity /
                                                      (h * h)));
}

/**
 * @brief Calculate the stabilization parameter for the transient Navier-Stokes
 * equations
 *
 * @param u_mag Magnitude of the velocity
 * @param kinematic_viscosity Kinematic viscosity
 * @param h Cell size; it should be calculated using the diameter of a sphere of
 * equal volume to that of the cell.
 * @param sdt Inverse of the time step (1/dt)
 *
 * @return Value of the stabilization parameter - tau
 */

inline double
calculate_navier_stokes_gls_tau_transient(const double u_mag,
                                          const double kinematic_viscosity,
                                          const double h,
                                          const double sdt)
{
  return 1. / std::sqrt(Utilities::fixed_power<2>(sdt) +
                        Utilities::fixed_power<2>(2. * u_mag / h) +
                        9 * Utilities::fixed_power<2>(4 * kinematic_viscosity /
                                                      (h * h)));
}

/**
 * @brief Compute the covariant metric tensor of an element from the inverse of
 * the mapping Jacobian.
 *
 * The covariant metric tensor is defined as
 * \f$ G_{ij} = \sum_k \frac{\partial \xi_k}{\partial x_i}
 *                       \frac{\partial \xi_k}{\partial x_j}
 *           = J^{-T} J^{-1} \f$
 * where \f$ \xi \f$ are the reference coordinates and \f$ x \f$ the real
 * coordinates. It carries the directional (anisotropic) size information of the
 * element and is used to build a metric-tensor-based stabilization parameter
 * tau (Bazilevs et al., variational multiscale methods) that remains robust for
 * deformed, stretched or high-aspect-ratio cells.
 *
 * @note deal.II's FEEvaluation::inverse_jacobian() and
 * FEValues::inverse_jacobian() return the inverse and transposed Jacobian
 * \f$ J^{-T} \f$: the \f$ (i,j) \f$ entry is \f$ \partial \xi_j / \partial x_i
 * \f$. Hence \f$ G = J^{-T} J^{-1} \f$ is obtained as \f$ K K^{T} \f$ (with
 * \f$ K \f$ the returned tensor), not \f$ K^{T} K \f$; the two differ for
 * sheared/curved cells.
 *
 * @note deal.II uses the reference cell \f$ [0,1]^{dim} \f$ whereas the VMS
 * literature (from which the constant @p c_i in
 * @ref calculate_navier_stokes_metric_tau is taken) uses \f$ [-1,1]^{dim} \f$.
 * The mapping derivatives therefore differ by a factor 2 per direction. In
 * addition, to reproduce the way Lethe accounts for the polynomial degree in
 * the isotropic definition (where the element size is divided by @p fe_degree),
 * the effective mapping derivatives are further multiplied by @p fe_degree. Both
 * effects are folded into the metric tensor, which is therefore scaled by
 * \f$ (2\,p)^2 \f$ with \f$ p = \f$ @p fe_degree. This keeps the definition
 * consistent and robust for higher-order elements.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Scalar type (e.g. double or VectorizedArray).
 *
 * @param inverse_jacobian Inverse and transposed mapping Jacobian
 * \f$ K_{ij} = \partial \xi_j / \partial x_i \f$, as returned by
 * FEEvaluation::inverse_jacobian() or FEValues::inverse_jacobian().
 * @param fe_degree Polynomial degree of the finite element.
 *
 * @return Covariant metric tensor \f$ G \f$.
 */
template <int dim, typename Number>
inline Tensor<2, dim, Number>
compute_covariant_metric_tensor(const Tensor<2, dim, Number> &inverse_jacobian,
                                const unsigned int            fe_degree)
{
  const double scaling = Utilities::fixed_power<2>(2. * fe_degree);
  return scaling * (inverse_jacobian * transpose(inverse_jacobian));
}

/**
 * @brief Compute the metric vector \f$ g_i = \sum_j \partial \xi_j / \partial
 * x_i \f$ used to define the grad-div (LSIC) stabilization parameter in the
 * metric-tensor formulation.
 *
 * As for @ref compute_covariant_metric_tensor, the vector is scaled by
 * \f$ 2\,p \f$ (with \f$ p = \f$ @p fe_degree) to account for the
 * \f$ [0,1]^{dim} \f$ vs \f$ [-1,1]^{dim} \f$ reference cell convention and for
 * the polynomial degree.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Scalar type (e.g. double or VectorizedArray).
 *
 * @param inverse_jacobian Inverse and transposed mapping Jacobian
 * \f$ K_{ij} = \partial \xi_j / \partial x_i \f$ (see
 * @ref compute_covariant_metric_tensor for the deal.II convention).
 * @param fe_degree Polynomial degree of the finite element.
 *
 * @return Metric vector \f$ g \f$.
 */
template <int dim, typename Number>
inline Tensor<1, dim, Number>
compute_metric_g_vector(const Tensor<2, dim, Number> &inverse_jacobian,
                        const unsigned int            fe_degree)
{
  Tensor<1, dim, Number> g;
  for (unsigned int i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < dim; ++j)
      g[i] += inverse_jacobian[i][j];
  return (2. * fe_degree) * g;
}

/**
 * @brief Calculate the metric-tensor-based (Bazilevs VMS) stabilization
 * parameter tau for the Navier-Stokes equations.
 *
 * \f$ \tau_M = \left( \mathrm{sdt}^2 + u \cdot G \cdot u + C_I \nu^2\, G:G
 * \right)^{-1/2} \f$
 *
 * where \f$ G \f$ is the covariant metric tensor
 * (@ref compute_covariant_metric_tensor), \f$ u \f$ the velocity, \f$ \nu \f$
 * the kinematic viscosity, and \f$ \mathrm{sdt} = 1/\Delta t \f$ (equal to 0 in
 * steady state). The constant \f$ C_I \f$ comes from an inverse estimate; the
 * value 36 corresponds to linear (Q1) elements.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Scalar type (e.g. double or VectorizedArray).
 *
 * @param velocity Advective velocity vector.
 * @param kinematic_viscosity Kinematic viscosity.
 * @param metric_tensor Covariant metric tensor \f$ G \f$.
 * @param sdt Inverse of the time step (1/dt); 0 for steady state.
 *
 * @return Value of the stabilization parameter tau.
 */
template <int dim, typename Number>
inline Number
calculate_navier_stokes_metric_tau(const Tensor<1, dim, Number> &velocity,
                                   const Number &kinematic_viscosity,
                                   const Tensor<2, dim, Number> &metric_tensor,
                                   const Number                 &sdt)
{
  constexpr double c_i = 36.;

  const Number u_G_u = velocity * (metric_tensor * velocity);
  const Number G_G   = scalar_product(metric_tensor, metric_tensor);

  return 1. / std::sqrt(Utilities::fixed_power<2>(sdt) + u_G_u +
                        c_i * Utilities::fixed_power<2>(kinematic_viscosity) *
                          G_G);
}

/**
 * @brief Calculate the metric-tensor-based grad-div (LSIC) stabilization
 * parameter, \f$ \tau_C = 1 / (\tau_M\, g \cdot g) \f$.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Scalar type (e.g. double or VectorizedArray).
 *
 * @param g_vector Metric vector \f$ g \f$ (@ref compute_metric_g_vector).
 * @param tau Momentum stabilization parameter \f$ \tau_M \f$
 * (@ref calculate_navier_stokes_metric_tau).
 *
 * @return Value of the grad-div (LSIC) stabilization parameter.
 */
template <int dim, typename Number>
inline Number
calculate_navier_stokes_metric_tau_lsic(const Tensor<1, dim, Number> &g_vector,
                                        const Number                 &tau)
{
  return 1. / (tau * (g_vector * g_vector));
}

/**
 * @brief Implements a Moe a posteriori shock capturing method to keep the field bounded. This limiter is only used when DG advection is selected for the Tracer.
 * The limiter is based on the implementation proposed by: Moe, Scott A.,
 * James A. Rossmanith, and David C. Seal. "A simple and effective high-order
 * shock-capturing limiter for discontinuous Galerkin methods." arXiv preprint
 * arXiv:1507.03024 (2015). https://doi.org/10.48550/arXiv.1507.03024
 *
 * The implementation follows the main idea of the paper, but we hardcode the
 * value of alpha in the article to 0. This is supposed to lead to a more
 * diffusive limiter, but in all honesty, all of our tries have shown that
 * this is already a not very dissipative limiter when using BDF time
 * integrator.
 *
 * @param dof_handler The DoFHandler used in the scalar simulation. This DOFHandler must have only a single component (a scalar equation).
 * @param locally_relevant_vector A solution vector that contains the locally_relevant solution. This vector is used to read the solution on the ghost cells and will be modified at the end.
 * @param locally_owned_vector A solution vector that contains the locally_owned solution. This vector will be used to write the new limited solution.
 */
template <int dim>
void
moe_scalar_limiter(const DoFHandler<dim> &dof_handler,
                   GlobalVectorType      &locally_relevant_vector,
                   GlobalVectorType      &locally_owned_vector);
#endif
