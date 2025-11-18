// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_solver_strategy_h
#define lethe_physics_solver_strategy_h

#include <core/parameters.h>

template <typename VectorType>
class PhysicsSolver;

/**
 * @brief Base class that works as an interface for all non-linear solvers for all non-linear systems of equations.
 *
 */
template <typename VectorType>
class PhysicsSolverStrategy
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] physics_solver A pointer to the physics solver to which the
   * non-linear solver is attached.
   *
   * @param[in] param Non-linear solver parameters as specified in the
   * simulation parameter file.
   *
   */
  PhysicsSolverStrategy(PhysicsSolver<VectorType>         *physics_solver,
                  const Parameters::NonLinearSolver &params);

  /**
   * @brief Destructor.
   *
   */
  virtual ~PhysicsSolverStrategy()
  {}

  /**
   * @brief Solve the non-linear system of equations.
   *
   */
  virtual void
  solve() = 0;


  /**
   * @brief Get the current newton iteration.
   *
   * @return Iteration number.
   *
   */
  inline unsigned int
  get_current_newton_iteration() const
  {
    return outer_iteration;
  }

protected:
  /**
   * @brief Physics solver for which we need a non-linear solver.
   *
   */
  PhysicsSolver<VectorType> *physics_solver;

  /**
   * @brief Non linear solver parameters.
   *
   */
  Parameters::NonLinearSolver params;

  /**
   * @brief Number of current Newton iteration.
   *
   */
  unsigned int outer_iteration;
};

template <typename VectorType>
PhysicsSolverStrategy<VectorType>::PhysicsSolverStrategy(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : physics_solver(physics_solver)
  , params(params)
  , outer_iteration(0)
{}

#endif
