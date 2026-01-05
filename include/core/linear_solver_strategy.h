// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_linear_solver_strategy_h
#define lethe_linear_solver_strategy_h

#include <core/physics_solver_strategy.h>

/**
 * @brief The linear solver strategy class is used to solve the physics as a linear problem.
 */
template <typename VectorType>
class LinearSolverStrategy : public PhysicsSolverStrategy<VectorType>
{
public:
  /**
   * @brief Constructor for the linear solver strategy.
   *
   * @param[in] physics_solver A pointer to the physics solver to which the
   * solving strategy is attached.
   *
   */
  LinearSolverStrategy(PhysicsSolver<VectorType> *physics_solver);


  /**
   * @brief Solve the linear system of equations.
   */
  void
  solve() override;
};

template <typename VectorType>
LinearSolverStrategy<VectorType>::LinearSolverStrategy(
  PhysicsSolver<VectorType> *physics_solver)
  : PhysicsSolverStrategy<VectorType>(physics_solver)
{}

template <typename VectorType>
void
LinearSolverStrategy<VectorType>::solve()
{
  PhysicsSolver<VectorType> *solver = this->physics_solver;

  solver->assemble_system_matrix();
  solver->setup_preconditioner();
  solver->assemble_system_rhs();
  solver->solve_linear_system();
}

#endif
