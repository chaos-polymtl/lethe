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

#include "parameters.h"
#include "physics_solver.h"

template <typename VectorType>
class NonLinearSolver
{
public:
  NonLinearSolver(PhysicsSolver<VectorType>* physics_solver,
                  const Parameters::NonLinearSolver &        params);

  void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                     time_stepping_method,
        const bool   is_initial_step,
        const double absolute_residual,
        const double relative_residual);

private:
  PhysicsSolver<VectorType>* physics_solver;
  Parameters::NonLinearSolver                params;
};

template <typename VectorType>
NonLinearSolver<VectorType>::NonLinearSolver(
  PhysicsSolver<VectorType>* physics_solver,
  const Parameters::NonLinearSolver &        params)
  : physics_solver(physics_solver)
  , params(params)
{}

template <typename VectorType>
void
NonLinearSolver<VectorType>::solve(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step,
  const double                                            absolute_residual,
  const double                                            relative_residual)
{
  double       current_res;
  double       last_res;
  bool         first_step      = is_initial_step;
  unsigned int outer_iteration = 0;
  last_res                     = 1.0;
  current_res                  = 1.0;
  while ((current_res > params.tolerance) &&
         outer_iteration < params.maxIterations)
    {
      physics_solver->set_evaluation_point(
        physics_solver->get_present_solution());

      physics_solver->assemble_matrix_rhs(time_stepping_method);

      if (outer_iteration == 0)
        {
          current_res = physics_solver->get_system_rhs().l2_norm();
          last_res    = current_res;
        }

      if (params.verbosity != Parameters::quiet)
        {
          physics_solver->get_ostream()
            << "Newton iteration: " << outer_iteration
            << "  - Residual:  " << current_res << std::endl;
        }

      /*
      physics_solver->solve_linear_system(first_step,
                          this->nsparam.linearSolver.minimum_residual,
                          this->nsparam.linearSolver.relative_residual);
                          */
      physics_solver->solve_linear_system(first_step,
                                          absolute_residual,
                                          relative_residual);

      for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
        {
          physics_solver->set_local_evaluation_point(
            physics_solver->get_present_solution());

          physics_solver->get_local_evaluation_point().add(
            alpha, physics_solver->get_newton_update());

          physics_solver->get_nonzero_constraints().distribute(
            physics_solver->get_local_evaluation_point());

          physics_solver->set_evaluation_point(
            physics_solver->get_local_evaluation_point());
            
          physics_solver->assemble_rhs(time_stepping_method);

          current_res = physics_solver->get_system_rhs().l2_norm();

          if (params.verbosity != Parameters::quiet)
            {
              physics_solver->get_ostream() << "\t\talpha = " << std::setw(6) << alpha
                          << std::setw(0) << " res = "
                          << std::setprecision(params.display_precision)
                          << current_res << std::endl;
            }

          if (current_res < 0.9 * last_res || last_res < params.tolerance)
            {
              break;
            }
        }


      physics_solver->set_present_solution(
        physics_solver->get_evaluation_point());
      last_res = current_res;
      ++outer_iteration;
    }
}

#endif
