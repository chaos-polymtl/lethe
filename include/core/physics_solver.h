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

#include <deal.II/lac/affine_constraints.h>

#include "non_linear_solver.h"
#include "parameters.h"

/**
 * An interface for all physics solver classes to derive from.
 */

template <typename VectorType>
class PhysicsSolver
{
public:
  PhysicsSolver(NonLinearSolver<VectorType> *non_linear_solver);

  ~PhysicsSolver()
  {
    delete non_linear_solver;
  }

  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) = 0;

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;

  virtual void
  solve_linear_system(const bool   initial_step,
                      const double absolute_residual,
                      const double relative_residual,
                      const bool   renewed_matrix = true) = 0;

  void
  solve_non_linear_system(
    const Parameters::SimulationControl::TimeSteppingMethod
               time_stepping_method,
    const bool first_iteration);

  virtual void
  apply_constraints()
  {
    nonzero_constraints.distribute(local_evaluation_point);
  }

  // Getters
  const Parameters::NonLinearSolver &
  get_params() const;
  const VectorType &
  get_system_rhs() const;
  const VectorType &
  get_evaluation_point() const;
  VectorType &
  get_local_evaluation_point();
  const VectorType &
  get_present_solution() const;
  const VectorType &
  get_newton_update() const;

  AffineConstraints<double> &
  get_nonzero_constraints();

  ConditionalOStream &
  get_ostream();

  // Setters
  void
  set_system_rhs(const VectorType &system_rhs);
  void
  set_evaluation_point(const VectorType &evaluation_point);
  void
  set_local_evaluation_point(const VectorType &local_evaluation_point);
  void
  set_present_solution(const VectorType &present_solution);
  void
  set_newton_update(const VectorType &newton_update);

protected:
  // TODO std::unique or std::shared pointer
  NonLinearSolver<VectorType> *non_linear_solver;

  VectorType system_rhs;
  VectorType evaluation_point;
  VectorType local_evaluation_point;
  VectorType present_solution;
  VectorType newton_update;

  AffineConstraints<double> nonzero_constraints;

  ConditionalOStream pcout;
};

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  NonLinearSolver<VectorType> *non_linear_solver)
  : non_linear_solver(non_linear_solver) // Default copy ctor
  , pcout({std::cout})
{}

template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_non_linear_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              first_iteration)
{
  this->non_linear_solver->solve(time_stepping_method, first_iteration);
}

template <typename VectorType>
const VectorType &
PhysicsSolver<VectorType>::get_system_rhs() const
{
  return system_rhs;
}

template <typename VectorType>
const VectorType &
PhysicsSolver<VectorType>::get_evaluation_point() const
{
  return evaluation_point;
}

template <typename VectorType>
VectorType &
PhysicsSolver<VectorType>::get_local_evaluation_point()
{
  return local_evaluation_point;
}

template <typename VectorType>
const VectorType &
PhysicsSolver<VectorType>::get_present_solution() const
{
  return present_solution;
}

template <typename VectorType>
const VectorType &
PhysicsSolver<VectorType>::get_newton_update() const
{
  return newton_update;
}

template <typename VectorType>
AffineConstraints<double> &
PhysicsSolver<VectorType>::get_nonzero_constraints()
{
  return nonzero_constraints;
}

template <typename VectorType>
ConditionalOStream &
PhysicsSolver<VectorType>::get_ostream()
{
  return pcout;
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::set_system_rhs(const VectorType &system_rhs)
{
  this->system_rhs = system_rhs;
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::set_evaluation_point(
  const VectorType &evaluation_point)
{
  this->evaluation_point = evaluation_point;
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::set_local_evaluation_point(
  const VectorType &local_evaluation_point)
{
  this->local_evaluation_point = local_evaluation_point;
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::set_present_solution(
  const VectorType &present_solution)
{
  this->present_solution = present_solution;
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::set_newton_update(const VectorType &newton_update)
{
  this->newton_update = newton_update;
}

#endif
