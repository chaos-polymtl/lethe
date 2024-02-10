/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef lethe_non_linear_solver_h
#define lethe_non_linear_solver_h

#include <core/parameters.h>

template <typename VectorType>
class PhysicsSolver;


/**
 * @brief Base class that works as an interface for all non-linear solvers for all non-linear systems of equations.
 *
 */
template <typename VectorType>
class NonLinearSolver
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
  NonLinearSolver(PhysicsSolver<VectorType>         *physics_solver,
                  const Parameters::NonLinearSolver &params);

  /**
   * @brief Destructor.
   *
   */
  virtual ~NonLinearSolver()
  {}

  /**
   * @brief Solve the non-linear system of equations.
   *
   * @param[in] is_initial_step Boolean variable that controls which constraints
   * are going to be applied to the equations depending on the time step.
   *
   */
  virtual void
  solve(const bool is_initial_step) = 0;


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
NonLinearSolver<VectorType>::NonLinearSolver(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : physics_solver(physics_solver)
  , params(params)
  , outer_iteration(0)
{}

#endif
