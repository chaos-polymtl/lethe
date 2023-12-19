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

*
 * Author: Simon Gauvin, Polytechnique Montreal, 2019
                           */

#ifndef lethe_non_linear_solver_h
#define lethe_non_linear_solver_h

#include <core/parameters.h>

template <typename VectorType>
class PhysicsSolver;


/**
 * @brief NonlinearSolver. Base class for all non-linear solver for non-linear systems of equations.
 * This class is an interface.
 */
template <typename VectorType>
class NonLinearSolver
{
public:
  /**
   * @brief Constructor for the NonLinearSolver.
   *
   * @param physics_solver A pointer to the physics solver to which the non-linear solver is attached
   *
   * @param param Non-linear solver parameters
   *
   */
  NonLinearSolver(PhysicsSolver<VectorType>         *physics_solver,
                  const Parameters::NonLinearSolver &params);

  virtual ~NonLinearSolver()
  {}

  /**
   * @brief Solve the non-linear system of equation.
   *
   * @param is_initial_step Boolean variable that controls which constraints are
   * going to be applied to the equations
   */
  virtual void
  solve(const bool is_initial_step) = 0;


  /**
   * @brief Return the current newton iteration.
   *
   */
  virtual unsigned int
  get_current_newton_iteration()
  {
    return outer_iteration;
  }

protected:
  PhysicsSolver<VectorType>  *physics_solver;
  Parameters::NonLinearSolver params;
  unsigned int                outer_iteration;
};

template <typename VectorType>
NonLinearSolver<VectorType>::NonLinearSolver(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : physics_solver(physics_solver)
  , params(params)
{}

#endif
