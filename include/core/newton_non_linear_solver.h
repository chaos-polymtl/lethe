/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Bruno Blais, Polytechnique Montreal, 2019 -
 */

#ifndef lethe_newton_non_linear_solver_h
#define lethe_newton_non_linear_solver_h

#include "non_linear_solver.h"

/**
 * @brief NewtonNonlinearSolver. Non-linear solver for non-linear systems of equations which uses a Newton
 * method with \alpha relaxation to ensure that the residual is monotonically
 * decreasing.
 */
template <typename VectorType>
class NewtonNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
  /**
   * @brief Constructor for the NewtonNonLinearSolver.
   *
   * @param physics_solver A pointer to the physics solver to which the non-linear solver is attached
   *
   * @param param Non-linear solver parameters
   *
   */
  NewtonNonLinearSolver(PhysicsSolver<VectorType> *        physics_solver,
                        const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equation.
   *
   * @param time_stepping_method Time stepping method being used. This is
   * required since the jacobian of the matrix is going to depend on the method
   * used
   *
   * @param is_initial_step Boolean variable that controls which constraints are
   * going to be applied to the equations
   *
   * @param force_matrix_renewal Boolean variable that controls if the Newton non-linear
   * solver will force the re-caculation of the jacobian matrix and the
   * reconstruction of the preconditioner at every iteration.
   */
  void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method,
        const bool is_initial_step,
        const bool force_matrix_renewal = true) override;
};

template <typename VectorType>
NewtonNonLinearSolver<VectorType>::NewtonNonLinearSolver(
  PhysicsSolver<VectorType> *        physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
NewtonNonLinearSolver<VectorType>::solve(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step,
  const bool)
{
  double       current_res;
  double       last_res;
  bool         first_step      = is_initial_step;
  unsigned int outer_iteration = 0;
  last_res                     = 1.0;
  current_res                  = 1.0;

  PhysicsSolver<VectorType> *solver = this->physics_solver;

  const int current_physics  = solver->get_current_physics();
  auto &    evaluation_point = solver->get_evaluation_point(current_physics);
  auto &    present_solution = solver->get_present_solution(current_physics);

  while ((current_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      evaluation_point = present_solution;

      solver->assemble_matrix_and_rhs(time_stepping_method);

      if (outer_iteration == 0)
        {
          auto &system_rhs = solver->get_system_rhs(current_physics);
          current_res      = system_rhs.l2_norm();
          last_res         = current_res;
        }

      if (this->params.verbosity != Parameters::Verbosity::quiet)
        {
          solver->pcout << "Newton iteration: " << outer_iteration
                        << "  - Residual:  " << current_res << std::endl;
        }

      solver->solve_linear_system(first_step);

      for (double alpha = 1.0; alpha > 1e-1; alpha *= 0.5)
        {
          auto &local_evaluation_point =
            solver->get_local_evaluation_point(current_physics);
          auto &newton_update    = solver->get_newton_update(current_physics);
          local_evaluation_point = present_solution;
          local_evaluation_point.add(alpha, newton_update);
          solver->apply_constraints();
          evaluation_point = local_evaluation_point;
          solver->assemble_rhs(time_stepping_method);

          auto &system_rhs = solver->get_system_rhs(current_physics);
          current_res      = system_rhs.l2_norm();

          if (this->params.verbosity != Parameters::Verbosity::quiet)
            {
              solver->pcout << "\t\talpha = " << std::setw(6) << alpha
                            << std::setw(0) << " res = "
                            << std::setprecision(this->params.display_precision)
                            << current_res << std::endl;
            }

          if (current_res < this->params.step_tolerance * last_res ||
              last_res < this->params.tolerance)
            {
              break;
            }
        }

      present_solution = evaluation_point;
      last_res         = current_res;
      ++outer_iteration;
    }
}

#endif
