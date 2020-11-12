#ifndef lethe_skip_newton_non_linear_solver
#define lethe_skip_newton_non_linear_solver

#include "non_linear_solver.h"

/**
 * @brief SkipNewtonNonlinearSolver. Non-linear solver for non-linear systems of equations which uses a Newton
 * method with \alpha relaxation to ensure that the residual is monotonically
 * decreasing. This non-linear solver only recalculates the jacobian matrix at a
 * given frequency.
 */
template <typename VectorType>
class SkipNewtonNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
  /**
   * @brief Constructor for the SkipNewtonNonLinearSolver.
   *
   * @param physics_solver A pointer to the physics solver to which the non-linear solver is attached
   *
   * @param param Non-linear solver parameters
   *
   * @param force_matrix_renewal A boolean variables that controls if the matrix and the preconditioner
   * will be forced to be recalculated even if the number of skipped iteration
   * has not been reached. This is generally used when the value of the time
   * step or the time stepping scheme changes.
   *
   */
  SkipNewtonNonLinearSolver(PhysicsSolver<VectorType> *        physics_solver,
                            const Parameters::NonLinearSolver &param);

  void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method,
        const bool is_initial_step,
        const bool force_matrix_renewal) override;


private:
  const Parameters::NonLinearSolver parameters;
  unsigned int                      consecutive_iters;
};

template <typename VectorType>
SkipNewtonNonLinearSolver<VectorType>::SkipNewtonNonLinearSolver(
  PhysicsSolver<VectorType> *        physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
  , parameters(params)
  , consecutive_iters(0)
{}

template <typename VectorType>
void
SkipNewtonNonLinearSolver<VectorType>::solve(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step,
  const bool                                              force_matrix_renewal)
{
  double       current_res;
  double       last_res;
  bool         first_step      = is_initial_step;
  unsigned int outer_iteration = 0;
  last_res                     = 1.0;
  current_res                  = 1.0;

  bool assembly_needed =
    consecutive_iters == 0 || is_initial_step || force_matrix_renewal;

  PhysicsSolver<VectorType> *solver     = this->physics_solver;
  auto &                     system_rhs = solver->get_system_rhs();

  while ((current_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      auto &evaluation_point = solver->get_evaluation_point();
      auto &present_solution = solver->get_present_solution();
      evaluation_point       = present_solution;

      if (assembly_needed)
        solver->assemble_matrix_and_rhs(time_stepping_method);

      else if (outer_iteration == 0)
        solver->assemble_rhs(time_stepping_method);

      if (outer_iteration == 0)
        {
          current_res = system_rhs.l2_norm();
          last_res    = current_res;
        }

      if (this->params.verbosity != Parameters::Verbosity::quiet)
        {
          solver->pcout << "Newton iteration: " << outer_iteration
                        << "  - Residual:  " << current_res << std::endl;
        }

      solver->solve_linear_system(first_step, assembly_needed);

      for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
        {
          auto &local_evaluation_point = solver->get_local_evaluation_point();
          auto &newton_update          = solver->get_newton_update();
          local_evaluation_point       = present_solution;
          local_evaluation_point.add(alpha, newton_update);
          solver->apply_constraints();
          evaluation_point = local_evaluation_point;
          solver->assemble_rhs(time_stepping_method);

          current_res = system_rhs.l2_norm();

          if (this->params.verbosity != Parameters::Verbosity::quiet)
            {
              solver->pcout << "\t\talpha = " << std::setw(6) << alpha
                            << std::setw(0) << " res = "
                            << std::setprecision(this->params.display_precision)
                            << current_res << std::endl;
            }

          if (current_res < 0.9 * last_res || last_res < this->params.tolerance)
            {
              break;
            }
        }

      present_solution = evaluation_point;
      last_res         = current_res;
      ++outer_iteration;
      assembly_needed = false;
    }
  if (!force_matrix_renewal)
    {
      consecutive_iters++;
      consecutive_iters = consecutive_iters % parameters.skip_iterations;
    }
}

#endif
