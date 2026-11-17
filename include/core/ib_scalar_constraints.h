// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_ib_scalar_constraints_h
#define lethe_ib_scalar_constraints_h

#include <core/ib_particle.h>
#include <core/parameters.h>
#include <core/vector.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <cstdint>
#include <vector>

using namespace dealii;

/**
 * @brief Impose a boundary condition on the surface of the particles of the
 * sharp-interface immersed boundary method, for a scalar field discretized
 * with continuous (FE_Q) elements.
 *
 * The method follows the sharp-interface immersed boundary method of the
 * fluid dynamics solver (L. Barbeau, S. Étienne, C. Béguin and B. Blais,
 * Computers & Fluids, 239, 105415, 2022): the equation of the degrees of
 * freedom located next to the surface of the particles is replaced by a
 * polynomial stencil (see IBStencil). The stencil extrapolates the solution of
 * a neighbouring fluid cell along the normal to the surface, where the
 * boundary condition is imposed. The boundary condition imposes either the
 * value of the field (dirichlet type), or the flux into the particles
 * \f$ D \partial u / \partial n \f$ (neumann and robin types), where
 * \f$ D \f$ is the diffusivity and \f$ n \f$ the normal pointing out of the
 * particle (see Parameters::IBTracerBoundaryCondition).
 *
 * Each locally owned or ghost cell is classified as:
 * - a fluid cell, if all its support points are outside the particles;
 * - a cell inside a particle, if all its support points are inside the
 *   particles;
 * - a cut cell otherwise. For a finite element of degree 1, a fluid cell whose
 *   degrees of freedom all belong to cut cells or are constrained is also
 *   treated as a cut cell, as in the fluid dynamics solver, since its equation
 *   would only involve constrained degrees of freedom.
 *
 * Each degree of freedom is then:
 * - a fluid degree of freedom if it only belongs to fluid cells. Its equation
 *   is assembled by the solver on the fluid cells.
 * - an interface degree of freedom if it belongs to at least one fluid cell
 *   and one cut cell. Its equation is replaced by the stencil that imposes
 *   the boundary condition.
 * - an extension degree of freedom otherwise, that is if it only belongs to
 *   cut cells and to cells inside the particles. Its value is a harmonic
 *   extension of the solution: a Laplace equation is assembled on the cut
 *   cells and on the cells inside the particles. The equations of the fluid
 *   and interface degrees of freedom do not involve the extension degrees of
 *   freedom, so the extension does not modify the solution in the fluid. It
 *   however gives a value consistent with the solution to the degrees of
 *   freedom that a moving particle uncovers.
 *
 * When no valid stencil is found for an interface degree of freedom, the
 * boundary condition is imposed to first order at the degree of freedom: the
 * dirichlet value is imposed on the degree of freedom, and the flux is
 * computed with the value and the gradient of the solution at the degree of
 * freedom, in a fluid cell that contains it.
 *
 * The solver must only assemble its equation on the fluid cells (see
 * is_fluid_cell()), then call assemble_matrix() and assemble_rhs().
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class IBScalarConstraints
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] p_stencil_degree Polynomial degree of the immersed boundary
   * stencil (see IBStencil).
   *
   * @param[in] p_length_ratio Length ratio of the immersed boundary stencil
   * (see IBStencil).
   *
   * @param[in] p_enable_extrapolation If false, the boundary condition is
   * imposed to first order on every interface degree of freedom, without
   * stencil. This is intended for debugging purposes.
   *
   * @param[in] p_boundary_condition Boundary condition imposed on the surface
   * of the particles. Its type must not be none. Its functions are shared with
   * the caller, which must set their time before assemble_rhs() is called.
   */
  IBScalarConstraints(
    const unsigned int                                p_stencil_degree,
    const double                                      p_length_ratio,
    const bool                                        p_enable_extrapolation,
    const Parameters::IBTracerBoundaryCondition<dim> &p_boundary_condition);

  /**
   * @brief Classify the cells and the degrees of freedom with respect to the
   * particles, and compute the stencil of every locally owned interface degree
   * of freedom. This function must be called after the degrees of freedom are
   * distributed and every time the particles move.
   *
   * @param[in] dof_handler DoFHandler of the scalar field. Its finite element
   * must be scalar, continuous and have support points.
   *
   * @param[in] geometry_dof_handler DoFHandler, built on the same
   * triangulation, whose cells are given as cell guesses to the shapes of the
   * particles. It is the DoFHandler on which the shapes compute their
   * precalculations, that is the DoFHandler of the fluid dynamics.
   *
   * @param[in] mapping Mapping of the scalar field.
   *
   * @param[in] particles Particles of the sharp-interface immersed boundary
   * method.
   *
   * @param[in] constraints Constraints of the Newton correction of the scalar
   * field (hanging nodes, Dirichlet and periodic boundary conditions).
   *
   * @param[in] p_diffusivity Diffusivity of the scalar field. It scales the
   * equations of the interface and extension degrees of freedom and, for the
   * neumann and robin types, it multiplies the normal derivative in the flux.
   * It must then be positive.
   *
   * @param[in] p_inverse_time_step Inverse of the time step, zero for a steady
   * simulation. It is only used to scale the equations of the interface and
   * extension degrees of freedom.
   */
  void
  update(const DoFHandler<dim>              &dof_handler,
         const DoFHandler<dim>              &geometry_dof_handler,
         const Mapping<dim>                 &mapping,
         const std::vector<IBParticle<dim>> &particles,
         const AffineConstraints<double>    &constraints,
         const double                        p_diffusivity,
         const double                        p_inverse_time_step);

  /**
   * @brief Return whether a cell is a fluid cell, on which the solver must
   * assemble its equation.
   *
   * @param[in] cell Locally owned or ghost cell of the DoFHandler given to
   * update().
   *
   * @return True if the cell is a fluid cell, false if it is cut by a particle
   * or inside a particle.
   */
  bool
  is_fluid_cell(
    const typename DoFHandler<dim>::active_cell_iterator &cell) const
  {
    AssertIndexRange(cell->active_cell_index(), cell_categories.size());
    return cell_categories[cell->active_cell_index()] == CellCategory::fluid;
  }

  /**
   * @brief Add the Laplace equation of the extension degrees of freedom to the
   * system matrix, then replace the equation of every interface degree of
   * freedom by its stencil. The solver must have assembled and compressed its
   * equation on the fluid cells before.
   *
   * @param[in] dof_handler DoFHandler given to update().
   *
   * @param[in] mapping Mapping given to update().
   *
   * @param[in] quadrature Quadrature used to assemble the Laplace equation.
   *
   * @param[in] constraints Constraints given to update().
   *
   * @param[in] evaluation_point Solution at which the Jacobian is evaluated.
   * It must have its ghost values. Only the robin type depends on it.
   *
   * @param[in,out] system_matrix System matrix (Jacobian) of the scalar field.
   */
  void
  assemble_matrix(const DoFHandler<dim>           &dof_handler,
                  const Mapping<dim>              &mapping,
                  const Quadrature<dim>           &quadrature,
                  const AffineConstraints<double> &constraints,
                  const GlobalVectorType          &evaluation_point,
                  TrilinosWrappers::SparseMatrix  &system_matrix) const;

  /**
   * @brief Add the residual of the Laplace equation of the extension degrees
   * of freedom to the right-hand side, then replace the right-hand side of
   * every interface degree of freedom by the residual of its stencil. The
   * right-hand side is the opposite of the residual, as for the Newton
   * correction. The solver must have assembled and compressed its equation on
   * the fluid cells before, and set the time of the functions of the boundary
   * condition.
   *
   * @param[in] dof_handler DoFHandler given to update().
   *
   * @param[in] mapping Mapping given to update().
   *
   * @param[in] quadrature Quadrature used to assemble the Laplace equation.
   *
   * @param[in] constraints Constraints given to update().
   *
   * @param[in] evaluation_point Solution at which the residual is evaluated.
   * It must have its ghost values.
   *
   * @param[in,out] system_rhs Right-hand side of the scalar field.
   */
  void
  assemble_rhs(const DoFHandler<dim>           &dof_handler,
               const Mapping<dim>              &mapping,
               const Quadrature<dim>           &quadrature,
               const AffineConstraints<double> &constraints,
               const GlobalVectorType          &evaluation_point,
               GlobalVectorType                &system_rhs) const;

  /**
   * @brief Return the flux of the scalar field into every particle, computed
   * from the residual of the equation that the solver assembles on the fluid
   * cells.
   *
   * The sum of the test functions of the interface degrees of freedom of a
   * particle is one on the boundary of the fluid region around this particle,
   * and zero on the rest of the boundary of the fluid region. The equation
   * tested with this sum is therefore a balance of the scalar over the fluid
   * cells next to the particle, whose only boundary term is the diffusive
   * flux into the particle (consistent flux). This flux is the sum of the
   * right-hand side of the interface degrees of freedom of the particle,
   * before their equation is replaced. It is the flux through the boundary of
   * the fluid region, which differs from the flux through the surface of the
   * particle by the accumulation, advection, reaction and sources in the cut
   * cells between them. These terms vanish in a steady diffusion problem
   * without sources, and their volume vanishes as the mesh is refined. A
   * positive flux enters the particle. The interface degrees of freedom that
   * are constrained, for example on the boundary of the domain, are not
   * counted, so the flux of a particle that touches the boundary of the domain
   * is underestimated.
   *
   * @param[in] fluid_cells_rhs Right-hand side of the equation assembled by the
   * solver on the fluid cells, at the solution, and compressed. The right-hand
   * side is the opposite of the residual, as for the Newton correction. The
   * equation of the interface degrees of freedom must not have been replaced
   * by assemble_rhs().
   *
   * @return Flux into every particle, in the order of the particles given to
   * update(), summed over all the processes.
   */
  std::vector<double>
  calculate_particle_fluxes(const GlobalVectorType &fluid_cells_rhs) const;

  /**
   * @brief Return the number of interface degrees of freedom whose equation is
   * replaced, summed over all the processes.
   *
   * @return Number of interface degrees of freedom.
   */
  types::global_dof_index
  get_n_interface_dofs() const
  {
    return n_interface_dofs;
  }

  /**
   * @brief Return the number of interface degrees of freedom on which the
   * boundary condition is imposed to first order because no valid stencil was
   * found, summed over all the processes.
   *
   * @return Number of interface degrees of freedom without stencil.
   */
  types::global_dof_index
  get_n_dofs_without_stencil() const
  {
    return n_dofs_without_stencil;
  }

private:
  /**
   * @brief Position of a cell with respect to the particles.
   */
  enum class CellCategory : std::uint8_t
  {
    /// All the support points of the cell are outside the particles.
    fluid,
    /// The cell is cut by the surface of a particle.
    cut,
    /// All the support points of the cell are inside the particles.
    inside
  };

  /**
   * @brief Role of a degree of freedom with respect to the particles.
   */
  enum class DoFCategory : std::uint8_t
  {
    /// The degree of freedom only belongs to fluid cells.
    fluid,
    /// The degree of freedom belongs to fluid and cut cells.
    interface,
    /// The degree of freedom only belongs to cut cells and inside cells.
    extension
  };

  /**
   * @brief Equation of a locally owned interface degree of freedom. The value
   * of the solution at the point of the surface closest to the degree of
   * freedom, \f$ \mathbf{x}_\Gamma \f$, is
   * \f$ u_\Gamma = \sum_m v_m u_m \f$, and its derivative along the normal
   * pointing out of the particle is
   * \f$ \partial_n u_\Gamma = \sum_m g_m u_m \f$, where \f$ v_m \f$ are the
   * value coefficients, \f$ g_m \f$ the normal derivative coefficients and
   * \f$ u_m \f$ the degrees of freedom of the stencil cell. The equation is
   * \f$ s \left( u_\Gamma - u_D(\mathbf{x}_\Gamma) \right) = 0 \f$ for the
   * dirichlet type and
   * \f$ s \left( J(u_\Gamma, \mathbf{x}_\Gamma) - D \partial_n u_\Gamma
   * \right) = 0 \f$ for the neumann and robin types, where \f$ s \f$ is the
   * scaling, \f$ u_D \f$ the imposed value and \f$ J \f$ the imposed flux.
   * With these signs, the diagonal coefficient of both equations is usually
   * positive, as for the equations assembled on the fluid cells.
   */
  struct InterfaceRow
  {
    /// Global index of the interface degree of freedom.
    types::global_dof_index dof_index;
    /// Index of the particle closest to the degree of freedom.
    unsigned int particle_index;
    /// Point of the surface of the particle closest to the degree of freedom.
    Point<dim> boundary_point;
    /// Global indices of the degrees of freedom involved in the equation.
    std::vector<types::global_dof_index> column_indices;
    /// Coefficients giving the value of the solution at the boundary point.
    std::vector<double> value_coefficients;
    /// Coefficients giving the derivative of the solution at the boundary
    /// point along the normal pointing out of the particle. It is empty for
    /// the dirichlet type.
    std::vector<double> normal_derivative_coefficients;
    /// Scaling of the equation, of the order of the diagonal of the matrix.
    double scaling;
  };

  /**
   * @brief Return the category of a locally relevant degree of freedom.
   *
   * @param[in] dof_index Global index of the degree of freedom.
   *
   * @return Category of the degree of freedom.
   */
  DoFCategory
  get_dof_category(const types::global_dof_index dof_index) const;

  /**
   * @brief Return the rate used to scale the equations of the degrees of
   * freedom of a cell: the sum of the inverse of the time step and of the
   * diffusivity divided by the square of the cell size. If both are zero, the
   * inverse of the square of the cell size is used.
   *
   * @param[in] cell Cell whose rate is computed.
   *
   * @return Rate of the cell, in units of inverse time.
   */
  double
  get_reference_rate(
    const typename DoFHandler<dim>::active_cell_iterator &cell) const;

  /// Polynomial degree of the immersed boundary stencil.
  const unsigned int stencil_degree;

  /// Length ratio of the immersed boundary stencil.
  const double length_ratio;

  /// If false, the boundary condition is imposed without stencil.
  const bool enable_extrapolation;

  /// Boundary condition imposed on the surface of the particles.
  const Parameters::IBTracerBoundaryCondition<dim> boundary_condition;

  /// True if the boundary condition imposes the flux (neumann and robin
  /// types), false if it imposes the value (dirichlet type).
  const bool imposes_flux;

  /// Diffusivity of the scalar field, used to scale the equations and in the
  /// flux.
  double diffusivity;

  /// Inverse of the time step (zero if steady), used to scale the equations.
  double inverse_time_step;

  /// MPI communicator of the DoFHandler given to update().
  MPI_Comm mpi_communicator;

  /// Number of particles given to update().
  unsigned int n_particles;

  /// Locally relevant degrees of freedom of the scalar field.
  IndexSet locally_relevant_dofs;

  /// Category of the cells, indexed by their active cell index. Only the
  /// entries of the locally owned and ghost cells are meaningful.
  std::vector<CellCategory> cell_categories;

  /// Category of the locally relevant degrees of freedom, indexed by their
  /// position in locally_relevant_dofs.
  std::vector<DoFCategory> dof_categories;

  /// Equations of the locally owned interface degrees of freedom that are not
  /// constrained.
  std::vector<InterfaceRow> interface_rows;

  /// Number of interface degrees of freedom, summed over all the processes.
  types::global_dof_index n_interface_dofs;

  /// Number of interface degrees of freedom without valid stencil, summed over
  /// all the processes.
  types::global_dof_index n_dofs_without_stencil;
};

#endif
