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

#ifndef LETHE_PHYSICSSOLVER
#define LETHE_PHYSICSSOLVER

#include "parameters.h"

/**
 * An interface for all physics solver classes to derive from.
 */

template <typename VectorType>
class PhysicsSolver
{
public:

  PhysicsSolver(const Parameters::NonLinearSolver& params);

  virtual void
  assemble_matrix_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) = 0;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;

  virtual void
  solve_linear_system(const bool       initial_step,
                      const double     absolute_residual,
                      const double     relative_residual) = 0;

  const Parameters::NonLinearSolver& get_params() const;
  VectorType& get_system_rhs() const;
  VectorType& get_evaluation_point() const;
  VectorType& get_local_evaluation_point() const;

protected:
  VectorType system_rhs;
  VectorType evaluation_point;
  VectorType local_evaluation_point;

  ConditionalOStream pcout;

private:
  const Parameters::NonLinearSolver& params;
};

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(const Parameters::NonLinearSolver& params)
  : pcout({std::cout})
  , params(params)
{}

template <typename VectorType>
const Parameters::NonLinearSolver& PhysicsSolver<VectorType>::get_params() const
{
  return params;
}

template <typename VectorType>
VectorType& PhysicsSolver<VectorType>::get_system_rhs() const
{
  return system_rhs;
}

template <typename VectorType>
VectorType& PhysicsSolver<VectorType>::get_evaluation_point() const
{
  return evaluation_point;
}

template <typename VectorType>
VectorType& PhysicsSolver<VectorType>::get_local_evaluation_point() const
{
  return local_evaluation_point;
}


#endif