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

#include "kinsol_newton_non_linear_solver.h"
#include "multiphysics.h"
#include "newton_non_linear_solver.h"
#include "non_linear_solver.h"
#include "parameters.h"

/**
 * This interface class is used to house all the common elements of physics
 * solver. A physics solver is an implementation of a linear or non-linear set
 * of physical equations. This interface is here to provide the families of
 * non-linear solvers with the necessary elements to allow for the solution of
 * the problems associated with a physics (block or not).
 */

template <typename VectorType>
class PhysicsSolver
{
public:
  /**
   * @brief PhysicsSolver
   * @param non_linear_solver_parameters A set of parameters that will be used to construct the non-linear solver
   * @param p_number_physic_total Indicates the number of physics solved
   * default value = 1, meaning only a single physics is solved
   */
  PhysicsSolver(const Parameters::NonLinearSolver non_linear_solver_parameters);

  virtual ~PhysicsSolver()
  {
    delete non_linear_solver;
  }

  /**
   * @brief Call for the assembly of the matrix
   */
  virtual void
  assemble_system_matrix(const Parameters::SimulationControl::TimeSteppingMethod
                           time_stepping_method) = 0;

  /**
   * @brief Call for the assembly of right-hand side
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_system_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) = 0;


  /**
   * @brief Call for the solution of the linear system of equation using a strategy appropriate
   * to the physics
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated or not
   */
  virtual void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) = 0;

  void
  solve_non_linear_system(
    const Parameters::SimulationControl::TimeSteppingMethod
               time_stepping_method,
    const bool first_iteration);

  virtual void
  apply_constraints()
  {
    auto &nonzero_constraints    = get_nonzero_constraints();
    auto &local_evaluation_point = get_local_evaluation_point();
    nonzero_constraints.distribute(local_evaluation_point);
  }


  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   * All methods derived from this physics must provide these elements. These
   * are the key ingredients which enable Lethe to solve a physics
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

  // attributes
  // TODO std::unique or std::shared pointer
  ConditionalOStream pcout;

private:
  NonLinearSolver<VectorType> *non_linear_solver;
};

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  const Parameters::NonLinearSolver non_linear_solver_parameters)
  : pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
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
      default:
        break;
    }
}

// solver method
template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_non_linear_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              first_iteration)
{
  // BB IMPORTANT
  // for (unsigned int iphys = 0; iphys < 1; iphys++)
  {
    this->non_linear_solver->solve(time_stepping_method, first_iteration);
  }
}
#endif
