// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_subequations_solver_h
#define lethe_physics_subequations_solver_h

#include <solvers/simulation_parameters.h>

#include <deal.II/numerics/data_out.h>

/**
 * @brief IDs associated to the different subequations solved in Lethe.
 */
enum SubequationsID : unsigned int
{
  /// VOF phase fraction gradient L2 projection
  phase_gradient_projection = 0
};

/**
 * @brief Base class to solve subequations solved inside a physics (or auxiliary
 * physics) that is not part of the main equations set. It contains all the
 * common elements of a subequation solver.
 */
class PhysicsSubequationsSolverBase
{
public:
  /**
   * @brief Set up the DofHandler and the degrees of freedom associated with
   * the equation to solve.
   */
  virtual void
  setup_dofs() = 0;

  /**
   * @brief Assemble the matrix.
   */
  virtual void
  assemble_system_matrix() = 0;

  /**
   * @brief Assemble the right-hand side (rhs).
   */
  virtual void
  assemble_system_rhs() = 0;

  /**
   * @brief Assemble and solve the system.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  virtual void
  solve(const bool &is_post_mesh_adaptation) = 0;
};

/**
 * @brief Linear subequations solved inside a physic (or auxiliary physic) that
 * is not part of the main set of equations.
 */

class PhysicsLinearSubequationsSolver : public PhysicsSubequationsSolverBase
{
public:
  /**
   * @brief Constructor for physics subequations that require the use of a
   * linear solver within the span of the physics resolution.
   *
   * @param[in] pcout Parallel cout used to print the information.
   */
  PhysicsLinearSubequationsSolver(const ConditionalOStream &pcout)
    : pcout(pcout)
  {}

  /**
   * @brief Solve the linear system associated with the equation to solve, when
   * the equation is already linear.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  virtual void
  solve_linear_system_and_update_solution(
    const bool &is_post_mesh_adaptation = false) = 0;

protected:
  ConditionalOStream pcout;
};

#endif
