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

#include <core/multiphysics.h>
#include <core/non_linear_solver.h>

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
  NewtonNonLinearSolver(PhysicsSolver<VectorType>         *physics_solver,
                        const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equation.
   *
   * @param is_initial_step Boolean variable that controls which constraints are
   * going to be applied to the equations
   */
  void
  solve(const bool is_initial_step) override;
};

template <typename VectorType>
NewtonNonLinearSolver<VectorType>::NewtonNonLinearSolver(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
NewtonNonLinearSolver<VectorType>::solve(const bool is_initial_step)
{
  double       global_res;
  double       current_res;
  double       global_res;
  double       last_res;
  bool         first_step      = is_initial_step;
  unsigned int outer_iteration = 0;
  last_res                     = 1e6;
  current_res                  = 1e6;
  global_res                   = 1e6;

  // current_res and global_res are different as one is defined based on the l2
  // norm of the residual vector (current_res) and the other (global_res) is
  // defined by the physical solver and may differ from the l2_norm of the
  // residual vector. Only the global_res is compared to the tolerance in order
  // to evaluate if the nonlinear system is solved. Only current_res is used for
  // the alpha scheme as this scheme only monitors the convergence of the
  // non-linear system of equation (the matrix problem).

  PhysicsSolver<VectorType> *solver = this->physics_solver;

  auto &evaluation_point = solver->get_evaluation_point();
  auto &present_solution = solver->get_present_solution();

  while ((global_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      evaluation_point = present_solution;

      solver->assemble_system_matrix();
      if (this->params.force_rhs_calculation || outer_iteration == 0)
        solver->assemble_system_rhs();

      if (outer_iteration == 0)
        {
          auto &system_rhs = solver->get_system_rhs();
          current_res      = system_rhs.l2_norm();
          last_res         = current_res;
        }

      if (this->params.verbosity != Parameters::Verbosity::quiet)
        {
          solver->pcout << "Newton iteration: " << outer_iteration
                        << "  - Residual:  " << current_res << std::endl;
        }

      solver->solve_linear_system(first_step);
      double last_alpha_res = current_res;

      unsigned int alpha_iter = 0;
      for (double alpha = 1.0; alpha > 1e-1; alpha *= 0.5)
        {
          auto &local_evaluation_point = solver->get_local_evaluation_point();
          auto &newton_update          = solver->get_newton_update();
          local_evaluation_point       = present_solution;
          local_evaluation_point.add(alpha, newton_update);
          solver->apply_constraints();
          evaluation_point = local_evaluation_point;
          solver->assemble_system_rhs();

          auto &system_rhs = solver->get_system_rhs();
          current_res      = system_rhs.l2_norm();

          if (this->params.verbosity != Parameters::Verbosity::quiet)
            {
              solver->pcout << "\t\talpha = " << std::setw(6) << alpha
                            << std::setw(0) << " res = "
                            << std::setprecision(this->params.display_precision)
                            << current_res << std::endl;
            }

          // If it's not the first iteration of alpha check if the residual is
          // smaller then the last alpha iteration. If it's not smaller we fall
          // back to the last alpha iteration.
          if (current_res > last_alpha_res and alpha_iter != 0)
            {
              alpha                  = 2 * alpha;
              local_evaluation_point = present_solution;
              local_evaluation_point.add(alpha, newton_update);
              solver->apply_constraints();
              evaluation_point = local_evaluation_point;

              if (this->params.verbosity != Parameters::Verbosity::quiet)
                {
                  solver->pcout
                    << "\t\talpha value was kept at alpha = " << alpha
                    << " since alpha = " << alpha / 2
                    << " increased the residual" << std::endl;
                }
              current_res = last_alpha_res;
              break;
            }
          if (current_res < this->params.step_tolerance * last_res ||
              last_res < this->params.tolerance)
            {
              break;
            }
          last_alpha_res = current_res;
          alpha_iter++;
        }

      global_res       = solver->get_current_residual();
      present_solution = evaluation_point;
      last_res         = current_res;
      ++outer_iteration;
    }
}

#endif
