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

#include "newton_non_linear_solver.h"
#include "non_linear_solver.h"
#include "parameters.h"
#include "skip_newton_non_linear_solver.h"

/**
 * This interface class is used to house all the common elements of physics
 * solver. A physics solver is an implementation of a linear or non-linear set
 * of physical equations. This interface is here to provide the families of
 * non-linear solvers with the necessary elements to allow for the solution of
 * the problems associated with a set of physics.
 */

template <typename VectorType>
class PhysicsSolver
{
public:
  /**
   * @brief PhysicsSolver - Constructor that takes an existing non-linear solver
   * @param non_linear_solver Non-linear solver that will be used to drive the physics solver
   */
  PhysicsSolver(NonLinearSolver<VectorType> *non_linear_solver);

  /**
   * @brief PhysicsSolver
   * @param non_linear_solver_parameters A set of parameters that will be used to construct the non-linear solver
   */
  PhysicsSolver(Parameters::NonLinearSolver non_linear_solver_parameters);

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
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) = 0;

  void
  solve_non_linear_system(
    const Parameters::SimulationControl::TimeSteppingMethod
               time_stepping_method,
    const bool first_iteration,
    const bool force_matrix_renewal);

  virtual void
  apply_constraints()
  {
    nonzero_constraints.distribute(local_evaluation_point);
  }

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
  , pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
{}



template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  Parameters::NonLinearSolver non_linear_solver_parameters)
  : pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
{
  switch (non_linear_solver_parameters.solver)
    {
      case Parameters::NonLinearSolver::SolverType::newton:
        non_linear_solver =
          new NewtonNonLinearSolver<VectorType>(this,
                                                non_linear_solver_parameters);
        break;
      case Parameters::NonLinearSolver::SolverType::skip_newton:
        non_linear_solver = new SkipNewtonNonLinearSolver<VectorType>(
          this, non_linear_solver_parameters);
        break;
      default:
        break;
    }
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_non_linear_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              first_iteration,
  const bool                                              force_matrix_renewal)
{
  this->non_linear_solver->solve(time_stepping_method,
                                 first_iteration,
                                 force_matrix_renewal);
}

#endif
