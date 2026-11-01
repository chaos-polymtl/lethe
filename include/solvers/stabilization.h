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
 * Reference-cell scaling: Bazilevs et al. define xi on the bi-unit cell
 * [-1,1]^d, whereas deal.II uses [0,1]^d, so dxi/dx is half as large.
 * Furthermore, Lethe's element size used by the other stabilizations is h/k
 * (see compute_cell_diameter). Passing reference_scaling = 2k multiplies dxi/dx
 * by that factor (G by (2k)^2 and g by 2k), so that tau_M -> h/(2|u|) and
 * tau_C -> |u|h/2 in the advective limit, consistently with the pspg_supg and
 * gls stabilizations.
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] inverse_jacobian J^{-T} as returned by inverse_jacobian(q), with
 *   entry [i][j] = dxi_j/dx_i.
 * @param[in] reference_scaling Scaling factor applied to dxi/dx (2k in Lethe,
 *   with k the FE degree of the velocity).
 * @param[out] metric_tensor The contravariant metric tensor G (eq. 66).
 * @param[out] metric_vector The metric vector g (eq. 69).
 */
template <int dim, typename Number>
inline void
compute_metric_tensor(const Tensor<2, dim, Number> &inverse_jacobian,
                      const double                  reference_scaling,
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

  // Map dxi/dx from the deal.II [0,1]^d cell to the bi-unit cell of Bazilevs
  // et al. and account for the FE degree (see reference-cell scaling above).
  metric_tensor *= reference_scaling * reference_scaling;
  metric_vector *= reference_scaling;
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
 * element-wise inverse estimate (Johnson; see paper text after eq. 70). Since
 * G already accounts for the FE degree (see compute_metric_tensor), Lethe uses
 * the constant C_I = 9, which recovers (in 1D) the viscous limit
 * tau_M = h^2/(12 nu) of the pspg_supg and gls stabilizations.
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
 * @brief Compute the stabilization parameter tau_t of the dynamic RBVMS
 * fine-scale velocity.
 *
 * The dynamic fine-scale velocity satisfies d(u')/dt + u'/tau_s = -r_M, where
 * tau_s is the steady stabilization parameter, i.e. calculate_rbvms_tau with
 * four_over_dt_squared = 0. Discretizing d(u')/dt with the same scheme as the
 * coarse-scale velocity, d(u')/dt ~ alpha u' + m', where m' gathers the
 * contributions of the previous fine-scale velocities, gives
 *
 *   u' = -tau_t (r_M + m'),   tau_t = ( alpha + 1/tau_s )^(-1)
 *
 * with alpha = beta_0 for BDF schemes and alpha = 0 for steady problems, in
 * which case tau_t = tau_s.
 *
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] tau_static The steady stabilization parameter tau_s.
 * @param[in] alpha The coefficient of the current fine-scale velocity in the
 *   time derivative (beta_0 for BDF schemes, 0 for steady problems).
 * @return The stabilization parameter tau_t.
 */
template <typename Number>
inline Number
calculate_rbvms_dynamic_tau(const Number tau_static, const double alpha)
{
  return 1.0 / (alpha + 1.0 / tau_static);
}

/**
 * @brief Assemble the RBVMS fine-scale cross-stress and Reynolds-stress
 * contributions to the velocity gradient residual (Bazilevs et al. 2007,
 * eq. 72), written in terms of the fine-scale velocity u'. These are the two
 * fine-scale terms that distinguish RBVMS from SUPG/PSPG. The returned tensor T
 * is meant to be added to the gradient residual that is later contracted with
 * the test-function gradient, i.e. it contributes
 * (T : ∇v) = sum_ij T[i][j] ∂v_i/∂x_j, so:
 *
 *   cross-stress    -( ∇v, u ⊗ u' )    -> -u_i u'_j
 *   Reynolds-stress -( ∇v, u' ⊗ u' )   -> -u'_i u'_j
 *
 * With the quasi-static fine-scale velocity u' = -τ_M r_M, these are the terms
 * +τ_M u_i r_M[j] and -τ_M^2 r_M[i] r_M[j] of eq. 72. Note the u factor sits
 * on the component (test) index i and u' on the derivative index j, which is
 * the transpose of the SUPG term -( u·∇v, u' ) = -u_j u'_i (u' on the
 * component index, u on the derivative index).
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] velocity The advecting velocity u.
 * @param[in] fine_scale_velocity The fine-scale velocity u'.
 * @return The gradient-residual contribution T[i][j].
 */
template <int dim, typename Number>
inline Tensor<2, dim, Number>
rbvms_fine_scale_terms(const Tensor<1, dim, Number> &velocity,
                       const Tensor<1, dim, Number> &fine_scale_velocity)
{
  Tensor<2, dim, Number> result;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      result[i][j] = -velocity[i] * fine_scale_velocity[j] -
                     fine_scale_velocity[i] * fine_scale_velocity[j];
  return result;
}

/**
 * @brief Linearization of rbvms_fine_scale_terms with respect to (δu, δu').
 * Applying the product rule to the cross-stress and Reynolds-stress terms of
 * eq. 72:
 *
 *   d(cross)    = -( δu_i u'_j + u_i δu'_j )
 *   d(Reynolds) = -( δu'_i u'_j + u'_i δu'_j )
 *
 * The fine-scale velocity is affine in the strong momentum residual r_M with a
 * frozen stabilization parameter (u' = -τ_M r_M for the quasi-static
 * subscale), so δu' = -τ_M δr_M, where δr_M is the linearization of r_M with
 * respect to the trial (δu, δp).
 *
 * @tparam dim Number of spatial dimensions.
 * @tparam Number Number type (e.g. double or VectorizedArray<double>).
 * @param[in] velocity Frozen advecting velocity u at the linearization point.
 * @param[in] delta_velocity Trial velocity increment δu.
 * @param[in] fine_scale_velocity Fine-scale velocity u' at the linearization
 *   point.
 * @param[in] delta_fine_scale_velocity Fine-scale velocity increment δu'.
 * @return The linearized gradient-residual contribution.
 */
template <int dim, typename Number>
inline Tensor<2, dim, Number>
rbvms_fine_scale_terms_linearization(
  const Tensor<1, dim, Number> &velocity,
  const Tensor<1, dim, Number> &delta_velocity,
  const Tensor<1, dim, Number> &fine_scale_velocity,
  const Tensor<1, dim, Number> &delta_fine_scale_velocity)
{
  Tensor<2, dim, Number> result;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      result[i][j] = -(delta_velocity[i] * fine_scale_velocity[j] +
                       velocity[i] * delta_fine_scale_velocity[j]) -
                     (delta_fine_scale_velocity[i] * fine_scale_velocity[j] +
                      fine_scale_velocity[i] * delta_fine_scale_velocity[j]);
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
