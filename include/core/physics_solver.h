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

#include "multiphysics.h"
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
   * @param p_number_physic_total Indicates the number of physics solved
   * default value = 1, meaning only a single physics is solved
   */
  // TODO delete the first constructor?
  PhysicsSolver(NonLinearSolver<VectorType> *non_linear_solver,
                unsigned int                 p_number_physic_total = 1);

  /**
   * @brief PhysicsSolver
   * @param non_linear_solver_parameters A set of parameters that will be used to construct the non-linear solver
   * @param p_number_physic_total Indicates the number of physics solved
   * default value = 1, meaning only a single physics is solved
   */
  PhysicsSolver(Parameters::NonLinearSolver non_linear_solver_parameters,
                unsigned int                p_number_physic_total = 1);

  ~PhysicsSolver()
  {
    delete non_linear_solver;
  }

  /**
   * @brief Call for the assembly of the matrix and the right-hand side
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) = 0;

  /**
   * @brief Call for the assembly of right-hand side
   *
   * @param time_stepping_method Time-Stepping method with which the assembly is called
   */
  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) = 0;


  /**
   * @brief Call for the solution of the linear system of equation
   *
   * @param initial_step Provides the linear solver with indication if this solution is the first
   * one for the system of equation or not
   *
   * @param renewed_matrix Indicates to the linear solve if the system matrix has been recalculated
   */
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
    auto &nonzero_constraints = get_nonzero_constraints(current_physics_id);
    auto &local_evaluation_point =
      get_local_evaluation_point(current_physics_id);


    nonzero_constraints.distribute(local_evaluation_point);
  }

  PhysicsID
  get_current_physics()
  {
    return current_physics_id;
  }

  /**
   * @brief Getter methods to get the private attributes for the physic currently solved
   *
   * @param number_physic_current Indicates the number associated with the physic currently solved
   * default value = 0, meaning only one physic, generally associated with the
   * flow equations, is solved by default
   */
  virtual VectorType &
  get_evaluation_point(const PhysicsID number_physic_current) = 0;
  virtual VectorType &
  get_local_evaluation_point(const PhysicsID number_physic_current) = 0;
  virtual VectorType &
  get_newton_update(const PhysicsID number_physic_current) = 0;
  virtual VectorType &
  get_present_solution(const PhysicsID number_physic_current) = 0;
  virtual VectorType &
  get_system_rhs(const PhysicsID number_physic_current) = 0;
  virtual AffineConstraints<double> &
  get_nonzero_constraints(const PhysicsID number_physic_current) = 0;

  // attributes
  // TODO std::unique or std::shared pointer
  ConditionalOStream pcout;

private:
  PhysicsID current_physics_id;

  NonLinearSolver<VectorType> *non_linear_solver;
};

// constructors
// TODO delete the first constructor?
template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  NonLinearSolver<VectorType> *non_linear_solver,
  unsigned int                 p_number_physic_total)
  : non_linear_solver(non_linear_solver) // Default copy ctor
  , pcout({std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0})
{}

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  Parameters::NonLinearSolver non_linear_solver_parameters,
  unsigned int                p_number_physic_total)
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

// solver method
template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_non_linear_system(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              first_iteration,
  const bool                                              force_matrix_renewal)
{
  // BB IMPORTANT
  // for (unsigned int iphys = 0; iphys < 1; iphys++)
  {
    this->current_physics_id = fluid_dynamics;
    this->non_linear_solver->solve(time_stepping_method,
                                   first_iteration,
                                   force_matrix_renewal);
  }
}
#endif
