// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_kinsol_newton_non_linear_solver_strategy_h
#define lethe_kinsol_newton_non_linear_solver_strategy_h

#include <core/parameters.h>
#include <core/physics_solver_strategy.h>

#include <deal.II/sundials/kinsol.h>

/**
 * @brief Non-linear Newton solver for non-linear systems of equations implemented
 * in the SUNDIALS suite, specifically the KINSOL package. This solver has
 * internal algorithms to determine the time step and decide whether to
 * reassemble the Jacobian matrix or not.
 */
template <typename VectorType>
class KinsolNewtonNonLinearSolverStrategy
  : public PhysicsSolverStrategy<VectorType>
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
  KinsolNewtonNonLinearSolverStrategy(PhysicsSolver<VectorType> *physics_solver,
                                      const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equations.
   */
  void
  solve() override;
};

template <typename VectorType>
KinsolNewtonNonLinearSolverStrategy<VectorType>::
  KinsolNewtonNonLinearSolverStrategy(PhysicsSolver<VectorType> *physics_solver,
                                      const Parameters::NonLinearSolver &params)
  : PhysicsSolverStrategy<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
KinsolNewtonNonLinearSolverStrategy<VectorType>::solve()
{
#ifdef DEAL_II_WITH_SUNDIALS

  PhysicsSolver<VectorType> *solver = this->physics_solver;

  auto &local_evaluation_point = solver->get_local_evaluation_point();

  VectorType &evaluation_point = solver->get_evaluation_point();
  VectorType &present_solution = solver->get_present_solution();

  typename SUNDIALS::KINSOL<VectorType>::AdditionalData additional_data;
  additional_data.function_tolerance            = this->params.tolerance;
  additional_data.maximum_non_linear_iterations = this->params.max_iterations;
  additional_data.step_tolerance                = this->params.tolerance;

  // Kinsol solution strategy, the default one is line search
  if (this->params.kinsol_strategy ==
      Parameters::NonLinearSolver::KinsolStrategy::normal_newton)
    {
      additional_data.strategy =
        SUNDIALS::KINSOL<VectorType>::AdditionalData::SolutionStrategy::newton;
    }
  else if (this->params.kinsol_strategy ==
           Parameters::NonLinearSolver::KinsolStrategy::line_search)
    {
      additional_data.strategy = SUNDIALS::KINSOL<
        VectorType>::AdditionalData::SolutionStrategy::linesearch;
    }
  else if (this->params.kinsol_strategy ==
           Parameters::NonLinearSolver::KinsolStrategy::picard)
    {
      additional_data.strategy =
        SUNDIALS::KINSOL<VectorType>::AdditionalData::SolutionStrategy::picard;
    }

  SUNDIALS::KINSOL<VectorType> nonlinear_solver(additional_data);

  nonlinear_solver.reinit_vector = [&](VectorType &x) {
    x.reinit(local_evaluation_point);
  };

  nonlinear_solver.residual = [&](const VectorType &evaluation_point_for_kinsol,
                                  VectorType       &residual) {
    solver->pcout << "Computing residual vector..." << std::endl;
    evaluation_point = evaluation_point_for_kinsol;
    solver->apply_constraints();
    solver->assemble_system_rhs();
    auto &current_residual = solver->get_system_rhs();
    residual               = current_residual;
    solver->pcout << "      -Residual: " << residual.l2_norm() << std::endl;
    return 0;
  };

  nonlinear_solver.setup_jacobian =
    [&](const VectorType &present_solution_for_kinsol,
        const VectorType & /*current_f*/) {
      solver->pcout << "Computing jacobian matrix..." << std::endl;
      evaluation_point = present_solution_for_kinsol;
      solver->assemble_system_matrix();
      solver->setup_preconditioner();
      return 0;
    };

  nonlinear_solver.solve_with_jacobian = [&](const VectorType & /* residual */,
                                             VectorType &dst,
                                             const double /* tolerance */) {
    solver->pcout << "Solving linear system..." << std::endl;
    solver->solve_linear_system();
    dst = solver->get_newton_update();
    return 0;
  };

  local_evaluation_point = present_solution;

  nonlinear_solver.solve(local_evaluation_point);
  present_solution = local_evaluation_point;
#else
  throw std::runtime_error(
    "Kinsol newton nonlinear solver requires DEAL_II to be compiled with SUNDIALS");
#endif // DEAL_II_WITH_SUNDIALS
}

#endif
