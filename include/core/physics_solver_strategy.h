// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_solver_strategy_h
#define lethe_physics_solver_strategy_h

#include <core/parameters.h>

#include <deal.II/base/exceptions.h>

template <typename VectorType>
class PhysicsSolver;

/**
 * @brief Base class that works as an interface for all solver strategies (either non-linear or linear) for all systems of equations.
 *
 */
template <typename VectorType>
class PhysicsSolverStrategy
{
public:
  /**
   * @brief Constructor for the non-linear solver strategies.
   *
   * @param[in] physics_solver A pointer to the physics solver to which the
   * solving strategy is attached.
   *
   * @param[in] param Non-linear solver parameters as specified in the
   * simulation parameter file.
   *
   */
  PhysicsSolverStrategy(PhysicsSolver<VectorType>         *physics_solver,
                        const Parameters::NonLinearSolver &params);

  /**
   * @brief Constructor for the linear solver strategy.
   *
   * @param[in] physics_solver A pointer to the physics solver to which the
   * solving strategy is attached.
   *
   */
  PhysicsSolverStrategy(PhysicsSolver<VectorType> *physics_solver);

  /**
   * @brief Destructor.
   *
   */
  virtual ~PhysicsSolverStrategy()
  {}

  /**
   * @brief Solve the system of equations.
   *
   */
  virtual void
  solve() = 0;


  /**
   * @brief Get the current newton iteration.
   *
   * @return Iteration number.
   *
   */
  inline unsigned int
  get_current_newton_iteration() const
  {
    return outer_iteration;
  }

protected:
  /**
   * @brief Skip the inner linear solve when the current Newton evaluation
   * already produced a freshly assembled RHS whose residual is below the
   * nonlinear tolerance.
   *
   * This shortcut is only valid when the RHS was assembled for the current
   * evaluation point. Otherwise the residual may still correspond to an older
   * state and skipping the solve would incorrectly reuse stale information.
   *
   * Solvers must opt into this behavior explicitly through
   * allow_skip_linear_solve_when_residual_is_below_tolerance(). When the
   * shortcut is actually taken, force_rhs_calculation must also be enabled so
   * that the converged residual comes from a fresh RHS assembled at the
   * current Newton evaluation point.
   *
   * @param[in] rhs_was_assembled_this_iteration Whether the current Newton
   * iteration assembled the RHS at the current evaluation point.
   * @param[in] rescale_metric Residual rescaling factor used by the nonlinear
   * solver strategy.
   * @param[in,out] current_res Current matrix residual norm.
   * @param[in,out] last_res Last accepted matrix residual norm.
   * @param[in,out] global_res Current nonlinear residual reported by the
   * physics solver.
   * @param[in,out] outer_iteration Current Newton iteration counter.
   *
   * @return true if the linear solve was skipped and the Newton state was
   * advanced consistently, false otherwise.
   */
  bool
  skip_linear_solve_if_fresh_rhs_is_already_converged(
    const bool    rhs_was_assembled_this_iteration,
    const double  rescale_metric,
    double       &current_res,
    double       &last_res,
    double       &global_res,
    unsigned int &outer_iteration);

  /**
   * @brief Physics solver for which we need a non-linear solver.
   *
   */
  PhysicsSolver<VectorType> *physics_solver;

  /**
   * @brief Non linear solver parameters.
   *
   */
  Parameters::NonLinearSolver params;

  /**
   * @brief Number of current Newton iteration.
   *
   */
  unsigned int outer_iteration;
};

template <typename VectorType>
PhysicsSolverStrategy<VectorType>::PhysicsSolverStrategy(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : physics_solver(physics_solver)
  , params(params)
  , outer_iteration(0)
{}

template <typename VectorType>
PhysicsSolverStrategy<VectorType>::PhysicsSolverStrategy(
  PhysicsSolver<VectorType> *physics_solver)
  : physics_solver(physics_solver)
{}

template <typename VectorType>
bool
PhysicsSolverStrategy<VectorType>::
  skip_linear_solve_if_fresh_rhs_is_already_converged(
    const bool    rhs_was_assembled_this_iteration,
    const double  rescale_metric,
    double       &current_res,
    double       &last_res,
    double       &global_res,
    unsigned int &outer_iteration)
{
  PhysicsSolver<VectorType> *solver = this->physics_solver;

  if (!solver->allow_skip_linear_solve_when_residual_is_below_tolerance())
    return false;

  if (!rhs_was_assembled_this_iteration)
    return false;

  const double assembled_res =
    solver->get_system_rhs().l2_norm() / rescale_metric;
  if (assembled_res > this->params.tolerance)
    return false;

  AssertThrow(
    this->params.force_rhs_calculation,
    ExcMessage(
      "Skipping the linear solve based on the assembled RHS residual requires "
      "'force rhs calculation = true' so the residual is freshly assembled on "
      "every Newton iteration."));

  current_res                 = assembled_res;
  solver->get_newton_update() = 0;
  global_res                  = solver->get_current_residual() / rescale_metric;
  solver->get_present_solution() = solver->get_evaluation_point();
  last_res                       = current_res;
  ++outer_iteration;

  return true;
}

#endif
