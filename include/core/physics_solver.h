/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2021 by the Lethe authors
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
 * ---------------------------------------------------------------------*/

#ifndef LETHE_PHYSICSSOLVER
#define LETHE_PHYSICSSOLVER

#include <core/inexact_newton_non_linear_solver.h>
#include <core/kinsol_newton_non_linear_solver.h>
#include <core/multiphysics.h>
#include <core/newton_non_linear_solver.h>
#include <core/non_linear_solver.h>
#include <core/parameters.h>

#include <deal.II/lac/affine_constraints.h>

/**
 * @brief Class that has all the common elements of a physics solver. It creates the nonlinear solver
 * as specified by the user using the parameters file and provides all the
 * necessary elements needed by the solver to solve a physics problem.
 *
 * @param non_linear_solver_parameters A set of parameters that will be used to construct the non-linear solver
 */
template <typename VectorType>
class PhysicsSolver
{
public:
  PhysicsSolver(const Parameters::NonLinearSolver non_linear_solver_parameters);

  virtual ~PhysicsSolver()
  {
    delete non_linear_solver;
  }

  /**
   * @brief assemble_system_matrix Assembles the matrix
   */
  virtual void
  assemble_system_matrix() = 0;

  /**
   * @brief assemble_system_rhs Assembles the rhs
   */
  virtual void
  assemble_system_rhs() = 0;


  /**
   * @brief solve_linear_system Solves the linear system of equations
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated or not
   */
  virtual void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) = 0;

  /**
   * @brief solve_non_linear_system Solves the non linear system of equations
   *
   * @param first_iteration Indicates whether it is the first iteration of the non linear solver
   */
  void
  solve_non_linear_system(const bool first_iteration);

  /**
   * @brief Applies constraints to a local_evaluation_point
   */
  virtual void
  apply_constraints()
  {
    auto &nonzero_constraints    = get_nonzero_constraints();
    auto &local_evaluation_point = get_local_evaluation_point();
    nonzero_constraints.distribute(local_evaluation_point);
  }

  /**
   * @brief Getter methods that give access to the private attributes of the physics being solved.
   * These methods must be provided by the physics as they are the key to solve
   * a problem using Lethe.
   */
  virtual VectorType &
  get_evaluation_point() = 0;
  virtual VectorType &
  get_local_evaluation_point() = 0;
  virtual VectorType &
  get_newton_update() = 0;
  virtual VectorType &
  get_present_solution() = 0;
  virtual VectorType &
  get_system_rhs() = 0;
  virtual AffineConstraints<double> &
  get_nonzero_constraints() = 0;

  /**
   * @brief Default way to evaluate the residual for the nonlinear solver.
   * Some application may use more complex evaluation of the residual and
   * override this method.
   */
  virtual double
  get_current_residual()
  {
    auto &system_rhs = get_system_rhs();
    return system_rhs.l2_norm();
  }
  virtual unsigned int
  get_current_newton_iteration()
  {
    return non_linear_solver->get_current_newton_iteration();
  }

  ConditionalOStream                                pcout;
  Parameters::SimulationControl::TimeSteppingMethod time_stepping_method;

private:
  NonLinearSolver<VectorType> *non_linear_solver;
};

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  const Parameters::NonLinearSolver non_linear_solver_parameters)
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
{
  switch (non_linear_solver_parameters.solver)
    {
      case Parameters::NonLinearSolver::SolverType::newton:
        non_linear_solver =
          new NewtonNonLinearSolver<VectorType>(this,
                                                non_linear_solver_parameters);
        break;
      case Parameters::NonLinearSolver::SolverType::kinsol_newton:
        non_linear_solver = new KinsolNewtonNonLinearSolver<VectorType>(
          this, non_linear_solver_parameters);
        break;
      case Parameters::NonLinearSolver::SolverType::inexact_newton:
        non_linear_solver = new InexactNewtonNonLinearSolver<VectorType>(
          this, non_linear_solver_parameters);
        break;
      default:
        break;
    }
}

template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_non_linear_system(const bool first_iteration)
{
  {
    this->non_linear_solver->solve(first_iteration);
  }
}
#endif
