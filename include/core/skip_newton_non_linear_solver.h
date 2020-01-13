#ifndef LETHE_SKIPNEWTONNONLINEARSOLVER
#define LETHE_SKIPNEWTONNONLINEARSOLVER

#include "non_linear_solver.h"

template <typename VectorType>
class SkipNewtonNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
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

  PhysicsSolver<VectorType> *solver = this->physics_solver;

  while ((current_res > this->params.tolerance) &&
         outer_iteration < this->params.max_iterations)
    {
      solver->evaluation_point = solver->present_solution;

      if (assembly_needed)
        solver->assemble_matrix_and_rhs(time_stepping_method);

      else if (outer_iteration == 0)
        solver->assemble_rhs(time_stepping_method);

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

      solver->solve_linear_system(first_step, assembly_needed);

      for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
        {
          solver->local_evaluation_point = solver->present_solution;
          solver->local_evaluation_point.add(alpha, solver->newton_update);
          solver->apply_constraints();
          solver->evaluation_point = solver->local_evaluation_point;
          solver->assemble_rhs(time_stepping_method);

          current_res = solver->system_rhs.l2_norm();

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

      solver->present_solution = solver->evaluation_point;
      last_res                 = current_res;
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
