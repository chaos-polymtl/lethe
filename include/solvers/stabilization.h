// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_stabilization_h
#define lethe_stabilization_h

#include <core/vector.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

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
 * The solution is rescaled around the cell average, which is obtained by
 * integrating the solution over the cell. The limiter consequently conserves
 * the integral of the solution for every polynomial degree. The extrema of the
 * solution within a cell are sampled at the support points of the finite
 * element and at the quadrature points. The cells that share a vertex through
 * a periodic boundary are considered as neighbors.
 *
 * @param dof_handler The DoFHandler used in the scalar simulation. This DOFHandler must have only a single component (a scalar equation) and its finite element must have support points.
 * @param mapping The mapping used in the scalar simulation.
 * @param cell_quadrature The quadrature used to integrate the solution over the cells. It should be the quadrature used to assemble the scalar equation so that the limiter conserves the same discrete integral.
 * @param locally_relevant_vector A solution vector that contains the locally_relevant solution. This vector is used to read the solution on the ghost cells and will be modified at the end.
 * @param locally_owned_vector A solution vector that contains the locally_owned solution. This vector will be used to write the new limited solution.
 */
template <int dim>
void
moe_scalar_limiter(const DoFHandler<dim> &dof_handler,
                   const Mapping<dim>    &mapping,
                   const Quadrature<dim> &cell_quadrature,
                   GlobalVectorType      &locally_relevant_vector,
                   GlobalVectorType      &locally_owned_vector);

/**
 * @brief Implements the hierarchical vertex-based a posteriori limiter of Kuzmin to keep the field bounded. This limiter is only used when DG advection is selected for the Tracer.
 * The limiter is based on: Kuzmin, Dmitri. "A vertex-based hierarchical slope
 * limiter for p-adaptive discontinuous Galerkin methods." Journal of
 * Computational and Applied Mathematics 233.12 (2010): 3077-3085.
 * https://doi.org/10.1016/j.cam.2009.05.028
 *
 * The solution within a cell is decomposed into its mean, a linear part and a
 * remainder which contains all the terms of higher degree. The linear part is
 * multiplied by a factor \f$\alpha_1\f$ and the remainder by a factor
 * \f$\alpha_2\f$, both between zero and one, which preserves the mean:
 * - \f$\alpha_1\f$ is the largest factor for which the linear part, evaluated
 *   at every vertex of the cell, remains between the minimum and the maximum
 *   of the means of the cells that share this vertex;
 * - \f$\alpha_2\f$ is obtained with the same procedure applied to every
 *   component of the gradient, which is linearly reconstructed using the
 *   second derivatives.
 *
 * The limiter is hierarchical because it sets \f$\alpha_1 =
 * \max(\alpha_1,\alpha_2)\f$. A smooth extremum, at which the gradient varies
 * smoothly although the solution is a local extremum, is consequently not
 * limited. This requires a degree of at least two. For a degree of one, there
 * is no second derivative and a single factor, calculated from the values of
 * the solution at the vertices, multiplies the deviation from the mean. The
 * limiter has no adjustable parameter. The solution is not limited in a cell
 * of which all the vertices are located on the boundaries of the domain, which
 * is the case of every cell of a mesh that is a single cell thick.
 *
 * The implementation differs from the article on the following points:
 * - The gradient and the second derivatives are those of the L2 projection of
 *   the solution onto the polynomials of degree two, instead of the
 *   derivatives at the centroid. The two are identical for the P2 elements of
 *   the article, but the projection is significantly more robust for higher
 *   degrees since it accounts for the solution over the whole cell.
 * - For degrees above two, the article limits every order of derivative
 *   separately. Here, all the terms above the linear one are multiplied by
 *   \f$\alpha_2\f$.
 * - The bounds of the components of the gradient are widened by a small
 *   fraction (0.1%) of the largest norm of the gradient among the cells that
 *   share the vertex. Along a direction in which the solution does not vary,
 *   the gradient and its bounds are only made of noise, and comparing them
 *   without this relaxation removes the higher-order part of smooth solutions.
 * - The vertices located on a boundary of the domain are not used to calculate
 *   the factors. The means of the cells around such a vertex do not bound a
 *   solution that varies monotonically in the direction normal to the
 *   boundary, and using them would limit smooth solutions next to every
 *   boundary.
 *
 * @param dof_handler The DoFHandler used in the scalar simulation. This DOFHandler must have only a single component (a scalar equation) and its finite element must have support points.
 * @param mapping The mapping used in the scalar simulation.
 * @param cell_quadrature The quadrature used to integrate the solution over the cells. It should be the quadrature used to assemble the scalar equation so that the limiter conserves the same discrete integral.
 * @param locally_relevant_vector A solution vector that contains the locally_relevant solution. This vector is used to read the solution on the ghost cells and will be modified at the end.
 * @param locally_owned_vector A solution vector that contains the locally_owned solution. This vector will be used to write the new limited solution.
 */
template <int dim>
void
kuzmin_scalar_limiter(const DoFHandler<dim> &dof_handler,
                      const Mapping<dim>    &mapping,
                      const Quadrature<dim> &cell_quadrature,
                      GlobalVectorType      &locally_relevant_vector,
                      GlobalVectorType      &locally_owned_vector);
#endif
