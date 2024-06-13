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
 */

#ifndef lethe_inexact_newton_non_linear_solver_h
#define lethe_inexact_newton_non_linear_solver_h

#include <core/non_linear_solver.h>

#include <iomanip>

/**
 * @brief Non-linear solver for non-linear systems of equations which uses an inexact
 * Newton method with \f$\alpha\f$ relaxation to ensure that the residual is
 * monotonically decreasing. It allows to reuse the matrix and avoid its
 * reassembly for every Newton step.
 *
 */
template <typename VectorType>
class InexactNewtonNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
  /**
   * @brief Constructor.
   *
   * @param[in] physics_solver A pointer to the physics solver to which the
   * non-linear solver is attached.
   *
   * @param[in] param Non-linear solver parameters as specified in the
   * simulation parameter file.
   *
   */
  InexactNewtonNonLinearSolver(PhysicsSolver<VectorType> *physics_solver,
                               const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equations.
   *
   * @param[in] is_initial_step Boolean variable that controls which constraints
   * are going to be applied to the equations depending on the time step.
   *
   */
  void
  solve(const bool is_initial_step) override;

private:
  bool matrix_requires_assembly;
};

template <typename VectorType>
InexactNewtonNonLinearSolver<VectorType>::InexactNewtonNonLinearSolver(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
  , matrix_requires_assembly(true)
{}

template <typename VectorType>
void
InexactNewtonNonLinearSolver<VectorType>::solve(const bool is_initial_step)
{
  double       global_res;
  double       current_res;
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

  if (!this->params.reuse_matrix)
    matrix_requires_assembly = true;

  while ((global_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      evaluation_point = present_solution;

      if (matrix_requires_assembly)
        {
          solver->assemble_system_matrix();
          solver->setup_preconditioner();
        }

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
              solver->pcout << "\talpha = " << std::setw(6) << alpha
                            << std::setw(0) << " res = "
                            << std::setprecision(this->params.display_precision)
                            << std::setw(6) << current_res;

              solver->output_newton_update_norms(
                this->params.display_precision);
            }

          // If it's not the first iteration of alpha check if the residual is
          // smaller then the last alpha iteration. If it's not smaller we fall
          // back to the last alpha iteration, which is generally the first one.
          if (alpha_iter != 0 and current_res > last_alpha_res)
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
              // If the current residual has decreased  by less than
              // the last residual * matrix_tolerance or more than one alpha
              // iteration was required, the jacobian approximation has become
              // insufficiently performant and, consequently, it is preferable
              // to renew the jacobian matrix
              if (current_res > this->params.matrix_tolerance * last_res ||
                  alpha_iter > 0)
                matrix_requires_assembly = true;
              else
                matrix_requires_assembly = false;

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

  // If the non-linear solver has not converged abort simulation if
  // abort_at_convergence_failure=true
  if ((global_res > this->params.tolerance) &&
      outer_iteration >= this->params.max_iterations &&
      this->params.abort_at_convergence_failure)
    {
      throw(std::runtime_error(
        "Stopping simulation because the non-linear solver has failed to converge"));
    }
}

#endif
