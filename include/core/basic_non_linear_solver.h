#ifndef LETHE_BASICNONLINEARSOLVER
#define LETHE_BASICNONLINEARSOLVER

#include "non_linear_solver.h"

template <typename VectorType>
class BasicNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
    BasicNonLinearSolver(PhysicsSolver<VectorType> *        physics_solver,
                         const Parameters::NonLinearSolver &params,
                         const double                       absolute_residual,
                         const double                       relative_residual);

    void
    solve(const Parameters::SimulationControl::TimeSteppingMethod
                    time_stepping_method,
            const bool is_initial_step) override;
};

template <typename VectorType>
BasicNonLinearSolver<VectorType>::BasicNonLinearSolver(
  PhysicsSolver<VectorType> *        physics_solver,
  const Parameters::NonLinearSolver &params,
  const double                       absolute_residual,
  const double                       relative_residual)
  : NonLinearSolver<VectorType>(physics_solver, params, absolute_residual, relative_residual)
{}

template <typename VectorType>
void
BasicNonLinearSolver<VectorType>::solve(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step)
{
  double       current_res;
  double       last_res;
  bool         first_step      = is_initial_step;
  unsigned int outer_iteration = 0;
  last_res                     = 1.0;
  current_res                  = 1.0;
  while ((current_res > this->params.tolerance) &&
         outer_iteration < this->params.maxIterations)
    {
      this->physics_solver->set_evaluation_point(
        this->physics_solver->get_present_solution());

      this->physics_solver->assemble_matrix_rhs(time_stepping_method);

      if (outer_iteration == 0)
        {
          current_res = this->physics_solver->get_system_rhs().l2_norm();
          last_res    = current_res;
        }

      if (this->params.verbosity != Parameters::quiet)
        {
          this->physics_solver->get_ostream()
            << "Newton iteration: " << outer_iteration
            << "  - Residual:  " << current_res << std::endl;
        }

      this->physics_solver->solve_linear_system(first_step,
                                          this->absolute_residual,
                                          this->relative_residual);

      for (double alpha = 1.0; alpha > 1e-3; alpha *= 0.5)
        {
          this->physics_solver->set_local_evaluation_point(
            this->physics_solver->get_present_solution());

          this->physics_solver->get_local_evaluation_point().add(
            alpha, this->physics_solver->get_newton_update());

          this->physics_solver->get_nonzero_constraints().distribute(
            this->physics_solver->get_local_evaluation_point());

          this->physics_solver->set_evaluation_point(
            this->physics_solver->get_local_evaluation_point());
            
          this->physics_solver->assemble_rhs(time_stepping_method);

          current_res = this->physics_solver->get_system_rhs().l2_norm();

          if (this->params.verbosity != Parameters::quiet)
            {
              this->physics_solver->get_ostream() << "\t\talpha = " << std::setw(6) << alpha
                          << std::setw(0) << " res = "
                          << std::setprecision(this->params.display_precision)
                          << current_res << std::endl;
            }

          if (current_res < 0.9 * last_res || last_res < this->params.tolerance)
            {
              break;
            }
        }


      this->physics_solver->set_present_solution(
        this->physics_solver->get_evaluation_point());
      last_res = current_res;
      ++outer_iteration;
    }
}

#endif