// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_stabilization_h
#define lethe_stabilization_h

#include <core/vector.h>

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

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
 * @brief Compute the element metric tensor G and metric vector g used by the
 * residual-based variational multiscale (RBVMS) stabilization from the inverse
 * of the mapping Jacobian.
 *
 * Following Bazilevs, Calo, Cottrell, Hughes, Reali & Scovazzi, "Variational
 * multiscale residual-based turbulence modeling for large eddy simulation of
 * incompressible flows", CMAME 197 (2007) 173-201:
 *
 *   G_ij = sum_k (dxi_k/dx_i)(dxi_k/dx_j)   (eq. 66, symmetric)
 *   g_i  = sum_j (dxi_j/dx_i)               (eq. 69)
 *
 * IMPORTANT (deal.II storage convention, verified against 9.7):
 * FEEvaluation::inverse_jacobian(q) returns J^{-T}, whose (i,j) entry contains
 * dxi_j/dx_i (see FEEvaluationData::inverse_jacobian docs: "columns refer to
 * reference space coordinates and rows to real cell coordinates"). This is the
 * transpose of a naive [d][e]=dxi_d/dx_e reading. With dxi_k/dx_i =
 * inv_j[i][k], the Bazilevs formulae become ROW operations on the returned
 * tensor:
 *
 *   G_ij = sum_k inv_j[i][k] * inv_j[j][k]   (eq. 66)
 *   g_i  = sum_j inv_j[i][j]                 (eq. 69)
 *
 * G is symmetric, so the transpose ambiguity only manifests on sheared cells;
 * it is checked explicitly by the rbvms_metric_tensor unit test.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] inverse_jacobian J^{-T} as returned by inverse_jacobian(q), with
 *   entry [i][j] = dxi_j/dx_i.
 * @param[out] metric_tensor The contravariant metric tensor G (eq. 66).
 * @param[out] metric_vector The metric vector g (eq. 69).
 */
template <int dim, typename Number>
inline void
compute_metric_tensor(const Tensor<2, dim, Number> &inverse_jacobian,
                      Tensor<2, dim, Number>       &metric_tensor,
                      Tensor<1, dim, Number>       &metric_vector)
{
  for (int i = 0; i < dim; ++i)
    {
      // g_i = sum_j dxi_j/dx_i = sum_j inv_j[i][j] (eq. 69, row sum)
      metric_vector[i] = inverse_jacobian[i][0];
      for (int j = 1; j < dim; ++j)
        metric_vector[i] += inverse_jacobian[i][j];

      // G_ij = sum_k dxi_k/dx_i * dxi_k/dx_j = sum_k inv_j[i][k]*inv_j[j][k]
      // (eq. 66, dot product of rows i and j)
      for (int j = 0; j < dim; ++j)
        {
          metric_tensor[i][j] = inverse_jacobian[i][0] * inverse_jacobian[j][0];
          for (int k = 1; k < dim; ++k)
            metric_tensor[i][j] +=
              inverse_jacobian[i][k] * inverse_jacobian[j][k];
        }
    }
}

/**
 * @brief Compute the RBVMS stabilization parameters tau_M (momentum) and tau_C
 * (continuity) from the element metric tensor G and metric vector g.
 *
 * Bazilevs et al. 2007, with the matrix-valued tau collapsed to the two scalars
 * of eq. 60:
 *
 *   tau_M = ( 4/dt^2 + u.G.u + C_I * nu^2 * (G:G) )^(-1/2)   (eq. 64)
 *   tau_C = 1 / ( tau_M * (g.g) )                            (eq. 65)
 *
 * with
 *   u.G.u = sum_ij u_i G_ij u_j   (eq. 68)
 *   G:G   = sum_ij G_ij G_ij      (eq. 67)
 *   g.g   = sum_i  g_i  g_i       (eq. 70)
 *
 * Both SUPG and PSPG use tau_M (eq. 72). C_I is the positive constant from an
 * element-wise inverse estimate (Johnson; see paper text after eq. 70); Lethe
 * uses the order-dependent value C_I = 3*k^2 (k = FE degree).
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] metric_tensor The contravariant metric tensor G (eq. 66).
 * @param[in] metric_vector The metric vector g (eq. 69).
 * @param[in] velocity The (advective) velocity u at the quadrature point.
 * @param[in] kinematic_viscosity The kinematic viscosity nu.
 * @param[in] four_over_dt_squared The transient term 4/dt^2 (0 for steady).
 * @param[in] c_i The inverse-estimate constant C_I.
 * @param[out] tau_momentum The momentum stabilization parameter tau_M (eq. 64).
 * @param[out] tau_continuity The continuity stabilization parameter tau_C
 *   (eq. 65).
 */
template <int dim, typename Number>
inline void
calculate_rbvms_tau(const Tensor<2, dim, Number> &metric_tensor,
                    const Tensor<1, dim, Number> &metric_vector,
                    const Tensor<1, dim, Number> &velocity,
                    const Number                  kinematic_viscosity,
                    const Number                  four_over_dt_squared,
                    const double                  c_i,
                    Number                       &tau_momentum,
                    Number                       &tau_continuity)
{
  Number u_G_u = 0.0; // eq. 68
  Number G_G   = 0.0; // eq. 67
  Number g_g   = 0.0; // eq. 70
  for (int i = 0; i < dim; ++i)
    {
      g_g += metric_vector[i] * metric_vector[i];
      for (int j = 0; j < dim; ++j)
        {
          u_G_u += velocity[i] * metric_tensor[i][j] * velocity[j];
          G_G += metric_tensor[i][j] * metric_tensor[i][j];
        }
    }

  // tau_M = ( 4/dt^2 + u.G.u + C_I * nu^2 * G:G )^(-1/2)   (eq. 64)
  tau_momentum =
    1.0 / std::sqrt(four_over_dt_squared + u_G_u +
                    c_i * kinematic_viscosity * kinematic_viscosity * G_G);

  // tau_C = 1 / ( tau_M * g.g )   (eq. 65)
  tau_continuity = 1.0 / (tau_momentum * g_g);
}

/**
 * @brief Assemble the RBVMS fine-scale cross-stress and Reynolds-stress
 * contributions to the velocity gradient residual (Bazilevs et al. 2007,
 * eq. 72). These are the two fine-scale terms that distinguish RBVMS from
 * SUPG/PSPG. The returned tensor T is meant to be added to the gradient
 * residual that is later contracted with the test-function gradient, i.e. it
 * contributes (T : ∇v) = sum_ij T[i][j] ∂v_i/∂x_j, so:
 *
 *   cross-stress     ( u·(∇v)^T, τ_M r_M )        -> +τ_M u_i r_M[j]
 *   Reynolds-stress -( ∇v, τ_M r_M ⊗ τ_M r_M )    -> -τ_M^2 r_M[i] r_M[j]
 *
 * The cross-stress index pattern is the transpose of the SUPG one (which places
 * r_M on the component index and u on the derivative index).
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] velocity The advecting velocity u.
 * @param[in] momentum_residual The strong momentum residual r_M (eq. 62).
 * @param[in] tau_momentum The momentum stabilization parameter τ_M (eq. 64).
 * @return The gradient-residual contribution T[i][j].
 */
template <int dim, typename Number>
inline Tensor<2, dim, Number>
rbvms_fine_scale_terms(const Tensor<1, dim, Number> &velocity,
                       const Tensor<1, dim, Number> &momentum_residual,
                       const Number                  tau_momentum)
{
  Tensor<2, dim, Number> result;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      result[i][j] = tau_momentum * velocity[i] * momentum_residual[j] -
                     tau_momentum * tau_momentum * momentum_residual[i] *
                       momentum_residual[j];
  return result;
}

/**
 * @brief Consistent frozen-τ linearization of rbvms_fine_scale_terms with
 * respect to (δu, δr_M). τ_M is treated as constant (read from the precomputed
 * table), so only the residual factors are differentiated. Applying the product
 * rule to the cross-stress and Reynolds-stress terms of eq. 72:
 *
 *   d(cross)    =  τ_M   ( δu_i r_M[j] + u_i δr_M[j] )
 *   d(Reynolds) = -τ_M^2 ( δr_M[i] r_M[j] + r_M[i] δr_M[j] )
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] velocity Frozen advecting velocity u^n.
 * @param[in] delta_velocity Trial velocity increment δu.
 * @param[in] momentum_residual Frozen strong residual r_M^n (eq. 62).
 * @param[in] delta_momentum_residual Linearized residual δr_M.
 * @param[in] tau_momentum The (frozen) momentum stabilization parameter τ_M.
 * @return The linearized gradient-residual contribution.
 */
template <int dim, typename Number>
inline Tensor<2, dim, Number>
rbvms_fine_scale_terms_linearization(
  const Tensor<1, dim, Number> &velocity,
  const Tensor<1, dim, Number> &delta_velocity,
  const Tensor<1, dim, Number> &momentum_residual,
  const Tensor<1, dim, Number> &delta_momentum_residual,
  const Number                  tau_momentum)
{
  Tensor<2, dim, Number> result;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      result[i][j] = tau_momentum * (delta_velocity[i] * momentum_residual[j] +
                                     velocity[i] * delta_momentum_residual[j]) -
                     tau_momentum * tau_momentum *
                       (delta_momentum_residual[i] * momentum_residual[j] +
                        momentum_residual[i] * delta_momentum_residual[j]);
  return result;
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
