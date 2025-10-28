// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_subequations_solver_h
#define lethe_physics_subequations_solver_h

#include <core/physics_solver.h>



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
   */
  virtual void
  solve() = 0;
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
  virtual ~PhysicsLinearSubequationsSolver() = default;

protected:
  /**
   * @brief Solve the linear system associated with the equation to solve, when
   * the equation is already linear.
   */
  virtual void
  solve_linear_system_and_update_solution() = 0;


  const ConditionalOStream pcout;
};

/**
 * @brief Non-linear subequations solved inside a physics (or auxiliary physics)
 * that are not part of the main set of equations.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the problem is solved.
 *
 * @tparam VectorType The Vector type used for the solvers.
 */
template <int dim, typename VectorType>
class PhysicsNonlinearSubequationsSolver : public PhysicsSubequationsSolverBase,
                                           public PhysicsSolver<VectorType>
{
public:
  /**
   * @brief Constructor for physics subequations that require the use of a
   * non-linear solver within the span of the physics solving.
   *
   * @param[in] non_linear_solver_parameters Set of parameters used to construct
   * the non-linear solver.
   *
   * @param[in] pcout Parallel cout used to print the information.
   */
  PhysicsNonlinearSubequationsSolver(
    const Parameters::NonLinearSolver non_linear_solver_parameters,
    const ConditionalOStream         &pcout)
    : PhysicsSolver<VectorType>(non_linear_solver_parameters)
    , pcout(pcout)
  {}

  /**
   * @brief Default destructor.
   */
  virtual ~PhysicsNonlinearSubequationsSolver() = default;

  /**
   * @brief Setup preconditioner.
   *
   * @note Not used for the physics subequations, but needed for the compilation
   * of the non-linear solver.
   */
  void
  setup_preconditioner() override{};

protected:
  const ConditionalOStream pcout;
};

#endif
