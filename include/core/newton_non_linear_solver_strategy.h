// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_newton_non_linear_solver_strategy_h
#define lethe_newton_non_linear_solver_strategy_h

#include <core/physics_solver_strategy.h>

/**
 * @brief Non-linear solver for non-linear systems of equations which uses a Newton
 * method with \f$\alpha\f$ relaxation to ensure that the residual is
 * monotonically decreasing.
 */
template <typename VectorType>
class NewtonNonLinearSolverStrategy : public PhysicsSolverStrategy<VectorType>
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
  NewtonNonLinearSolverStrategy(PhysicsSolver<VectorType> *physics_solver,
                                const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equations.
   *
   */
  void
  solve() override;
};

template <typename VectorType>
NewtonNonLinearSolverStrategy<VectorType>::NewtonNonLinearSolverStrategy(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : PhysicsSolverStrategy<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
NewtonNonLinearSolverStrategy<VectorType>::solve()
{
  double global_res;
  double current_res;
  double last_res;
  this->outer_iteration = 0;
  last_res              = 1e6;
  current_res           = 1e6;
  global_res            = 1e6;

  // current_res and global_res are different as one is defined based on the l2
  // norm of the residual vector (current_res) and the other (global_res) is
  // defined by the physical solver and may differ from the l2_norm of the
  // residual vector. Only the global_res is compared to the tolerance in order
  // to evaluate if the nonlinear system is solved. Only current_res is used for
  // the alpha scheme as this scheme only monitors the convergence of the
  // non-linear system of equation (the matrix problem).

  PhysicsSolver<VectorType> *solver = this->physics_solver;

  const double rescale_metric = solver->get_residual_rescale_metric();

  auto &evaluation_point = solver->get_evaluation_point();
  auto &present_solution = solver->get_present_solution();

  while ((global_res > this->params.tolerance) &&
         this->outer_iteration < this->params.max_iterations)
    {
      evaluation_point = present_solution;

      solver->assemble_system_matrix();

      if (!this->params.reuse_preconditioner || this->outer_iteration == 0)
        solver->setup_preconditioner();

      if (this->params.force_rhs_calculation || this->outer_iteration == 0)
        solver->assemble_system_rhs();

      if (this->outer_iteration == 0)
        {
          auto &system_rhs = solver->get_system_rhs();
          current_res      = system_rhs.l2_norm() / rescale_metric;
          last_res         = current_res;
        }

      if (this->params.verbosity != Parameters::Verbosity::quiet)
        {
          solver->pcout << "Newton iteration: " << this->outer_iteration
                        << "  - Residual:  " << current_res << std::endl;
        }

      solver->solve_linear_system();
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
          current_res      = system_rhs.l2_norm() / rescale_metric;

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
          // smaller than the last alpha iteration. If it's not smaller, we fall
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

      global_res       = solver->get_current_residual() / rescale_metric;
      present_solution = evaluation_point;
      last_res         = current_res;
      ++this->outer_iteration;
    }

  // If the non-linear solver has not converged abort simulation if
  // abort_at_convergence_failure=true
  if ((global_res > this->params.tolerance) &&
      this->outer_iteration >= this->params.max_iterations &&
      this->params.abort_at_convergence_failure)
    {
      throw(std::runtime_error(
        "Stopping simulation because the non-linear solver has failed to converge"));
    }
}

#endif
