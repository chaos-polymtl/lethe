// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_subequations_solver_h
#define lethe_physics_subequations_solver_h

#include <core/physics_solver.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/numerics/data_out.h>

enum SubequationsID : unsigned int
{
  phase_gradient_projection = 0
};

/**
 * @brief Sub-equations solved inside a physic (or auxiliary physic) that is not
 * part of the main equations set are solved through the
 * PhysicsLinearSubequationsSolver object. It contains all the common elements
 * of a physic solver. It creates the non-linear solver as specified by the user
 * using the parameters file and provides all the necessary elements needed by
 * the solver to solve a physics problem.
 *
 * @tparam dim Number of dimensions of the problem
 *
 * @tparam VectorType Type of vector container used to store solutions
 */

template <int dim, typename VectorType>
class PhysicsLinearSubequationsSolver
{
public:
  /**
   * @brief Constructor for physics sub-equations that require the use of a
   * linear and/or non-linear solver within the span of the physics resolution.
   *
   * @param[in] pcout Parallel cout used to print the information.
   */
  PhysicsLinearSubequationsSolver(
    const ConditionalOStream &pcout
    //    const Parameters::NonLinearSolver non_linear_solver_parameters
    )
    : //    : PhysicsSolver<VectorType>(non_linear_solver_parameters),
    pcout(pcout)
  {}

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
   * @brief Solve the linear system associated with the equation to solve, when
   * the equation is already linear.
   */
  virtual void
  solve_linear_system_and_update_solution() = 0;

  /**
   * @brief Assemble and solve linear system when the equation to solve is
   * linear without using the non-linear solver interface.
   */
  virtual void
  solve() = 0;

protected:
  ConditionalOStream pcout;
};

#endif
