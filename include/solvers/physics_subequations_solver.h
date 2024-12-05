// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_subequations_solver_h
#define lethe_physics_subequations_solver_h

#include <core/physics_solver.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/numerics/data_out.h>

/**
 * @brief Base class to solve subequations solved inside a physics (or auxiliary
 * physics) that are not part of the main equations set. It contains all the
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
   * @brief Assemble and solve the system.
   *
   * @param[in] is_post_mesh_adaptation Indicates if the equation is being
   * solved during post_mesh_adaptation(), for verbosity.
   */
  virtual void
  solve(const bool &is_post_mesh_adaptation) = 0;
};

/**
 * @brief Linear subequations solved inside a physics (or auxiliary physics)
 * that are not part of the main set of equations.
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
   * @brief Default destructor.
   */
  ~PhysicsLinearSubequationsSolver() = default;

protected:
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


  const ConditionalOStream pcout;
};

/**
 * @brief TODO AMISHGA
 * @tparam dim
 * @tparam VectorType
 */
template <int dim, typename VectorType>
class PhysicsNonlinearSubequationsSolver : public PhysicsSubequationsSolverBase,
                                           public PhysicsSolver<VectorType>
{
public:
  /**
   * @brief TODO AMISHGA
   *
   * @param non_linear_solver_parameters
   * @param pcout
   */
  PhysicsNonlinearSubequationsSolver(
    const Parameters::NonLinearSolver non_linear_solver_parameters,
    const ConditionalOStream         &pcout)
    : PhysicsSolver<VectorType>(non_linear_solver_parameters)
    , pcout(pcout)
  {}

  /**
   * TODO AMISHGA
   * @brief Set up preconditioner. Not used for the physics subequations, but
   * needed for the compilation of the non-linear solver.
   */
  void
  setup_preconditioner() override{};

protected:
  const ConditionalOStream pcout;
};

#endif
