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

#ifndef LETHE_NONLINEARSOLVER
#define LETHE_NONLINEARSOLVER

#include <deal.II/lac/affine_constraints.h>

#include "parameters.h"

template <typename VectorType>
class PhysicsSolver;

template <typename VectorType>
class NonLinearSolver
{
public:
  NonLinearSolver(PhysicsSolver<VectorType> *        physics_solver,
                  const Parameters::NonLinearSolver &params,
                  const double                       absolute_residual,
                  const double                       relative_residual);
  
  virtual ~NonLinearSolver() {}

  virtual void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method,
        const bool is_initial_step) = 0;

protected:
  PhysicsSolver<VectorType>* physics_solver;
  Parameters::NonLinearSolver                params;

  const double absolute_residual;
  const double relative_residual;
};

template <typename VectorType>
NonLinearSolver<VectorType>::NonLinearSolver(
  PhysicsSolver<VectorType> *        physics_solver,
  const Parameters::NonLinearSolver &params,
  const double                       absolute_residual,
  const double                       relative_residual)
  : physics_solver(physics_solver)
  , params(params)
  , absolute_residual(absolute_residual)
  , relative_residual(relative_residual)
{}

#endif
