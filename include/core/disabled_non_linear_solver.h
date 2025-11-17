// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_disabled_non_linear_solver_h
#define lethe_disabled_non_linear_solver_h

#include <core/non_linear_solver.h>

/**
 * @brief The linear physics class is used to solve the physics as a linear problem.
 */
template <typename VectorType>
class DisabledNonLinearSolver : public NonLinearSolver<VectorType>
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
  DisabledNonLinearSolver(PhysicsSolver<VectorType>         *physics_solver,
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
};

template <typename VectorType>
DisabledNonLinearSolver<VectorType>::DisabledNonLinearSolver(
  PhysicsSolver<VectorType>         *physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
DisabledNonLinearSolver<VectorType>::solve(const bool is_initial_step)
{

  PhysicsSolver<VectorType> *solver = this->physics_solver;


  if (this->params.verbosity != Parameters::Verbosity::quiet)
    {
      solver->pcout << "Linear physics, solving physics as a linear problem." 
                    << std::endl;
    }

  solver->assemble_system_matrix();
  solver->setup_preconditioner();
  solver->assemble_system_rhs();

  // The linear system is solved only once, therefore initial_step is always
  // true
  solver->solve_linear_system(true);
}

#endif
