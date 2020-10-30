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

  while ((current_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      solver->evaluation_point = solver->present_solution;

      solver->assemble_matrix_and_rhs(time_stepping_method);

      if (outer_iteration == 0)
        {
          current_res = solver->system_rhs.l2_norm();
          last_res    = current_res;
        }

      if (this->params.verbosity != Parameters::Verbosity::quiet)
        {
          solver->pcout << "Newton iteration: " << outer_iteration
                        << "  - Residual:  " << current_res << std::endl;
        }

      solver->solve_linear_system(first_step);

      for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
        {
          auto &local_evaluation_point = solver->get_local_evaluation_point();
          local_evaluation_point       = solver->present_solution;
          local_evaluation_point.add(alpha, solver->newton_update);
          solver->apply_constraints();
          solver->evaluation_point = local_evaluation_point;
          solver->assemble_rhs(time_stepping_method);

          current_res = solver->system_rhs.l2_norm();

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

      solver->present_solution = solver->evaluation_point;
      last_res                 = current_res;
      ++outer_iteration;
    }
}

#endif
