/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2021 -  by the Lethe authors
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

#ifndef kinsol_non_linear_solver_h
#define kinsol_non_linear_solver_h

#include <core/multiphysics.h>
#include <core/non_linear_solver.h>

#include <deal.II/sundials/kinsol.h>

/**
 * @brief KinsolNonlinearSolver. Interface to the non-linear newton solver for non-linear systems
 * of equations implemented in the SUNDIALS suite, specifically the KINSOL
 * package. This solver has internal algorithms to determine the time step and
 * decide whether to reassemble the Jacobian matrix or not.
 */
template <typename VectorType>
class KinsolNewtonNonLinearSolver : public NonLinearSolver<VectorType>
{
public:
  /**
   * @brief Constructor for the KinsolNewtonNonLinearSolver.
   *
   * @param physics_solver A pointer to the physics solver to which the non-linear solver is attached
   *
   * @param param Non-linear solver parameters
   *
   */
  KinsolNewtonNonLinearSolver(PhysicsSolver<VectorType> *        physics_solver,
                              const Parameters::NonLinearSolver &param);


  /**
   * @brief Solve the non-linear system of equation.
   *
   * @param time_stepping_method Time stepping method being used. This is
   * required since the jacobian of the matrix is going to depend on the method
   * used
   *
   * @param is_initial_step Boolean variable that controls which constraints are
   * going to be applied to the equations
   */
  void
  solve(const Parameters::SimulationControl::TimeSteppingMethod
                   time_stepping_method,
        const bool is_initial_step) override;
};

template <typename VectorType>
KinsolNewtonNonLinearSolver<VectorType>::KinsolNewtonNonLinearSolver(
  PhysicsSolver<VectorType> *        physics_solver,
  const Parameters::NonLinearSolver &params)
  : NonLinearSolver<VectorType>(physics_solver, params)
{}

template <typename VectorType>
void
KinsolNewtonNonLinearSolver<VectorType>::solve(
  const Parameters::SimulationControl::TimeSteppingMethod time_stepping_method,
  const bool                                              is_initial_step)
{
  (void)time_stepping_method;
  (void)is_initial_step;

  // double       current_res;
  // double       last_res;
  // bool         first_step      = is_initial_step;
  // unsigned int outer_iteration = 0;
  // last_res                     = 1e6;
  // current_res                  = 1e6;

  // PhysicsSolver<VectorType> *solver = this->physics_solver;

  // auto &evaluation_point = solver->get_evaluation_point();
  // auto &present_solution = solver->get_present_solution();

  // while ((current_res > this->params.tolerance) &&
  //        outer_iteration < this->params.max_iterations)
  //   {
  //     evaluation_point = present_solution;

  //     solver->assemble_matrix_and_rhs(time_stepping_method);

  //     if (outer_iteration == 0)
  //       {
  //         auto &system_rhs = solver->get_system_rhs();
  //         current_res      = system_rhs.l2_norm();
  //         last_res         = current_res;
  //       }

  //     if (this->params.verbosity != Parameters::Verbosity::quiet)
  //       {
  //         solver->pcout << "Newton iteration: " << outer_iteration
  //                       << "  - Residual:  " << current_res << std::endl;
  //       }

  //     solver->solve_linear_system(first_step);
  //     double last_alpha_res = current_res;

  //     for (double alpha = 1.0; alpha > 1e-1; alpha *= 0.5)
  //       {
  //         auto &local_evaluation_point =
  //         solver->get_local_evaluation_point(); auto &newton_update =
  //         solver->get_newton_update(); local_evaluation_point       =
  //         present_solution; local_evaluation_point.add(alpha, newton_update);
  //         solver->apply_constraints();
  //         evaluation_point = local_evaluation_point;
  //         solver->assemble_rhs(time_stepping_method);

  //         auto &system_rhs = solver->get_system_rhs();
  //         current_res      = system_rhs.l2_norm();

  //         if (this->params.verbosity != Parameters::Verbosity::quiet)
  //           {
  //             solver->pcout << "\t\talpha = " << std::setw(6) << alpha
  //                           << std::setw(0) << " res = "
  //                           <<
  //                           std::setprecision(this->params.display_precision)
  //                           << current_res << std::endl;
  //           }

  //         // If it's not the first iteration of alpha check if the residual
  //         is
  //         // smaller then the last alpha iteration. If it's not smaller we
  //         fall
  //         // back to the last alpha iteration.
  //         if (current_res > last_alpha_res and alpha < 0.99)
  //           {
  //             alpha                  = 2 * alpha;
  //             local_evaluation_point = present_solution;
  //             local_evaluation_point.add(alpha, newton_update);
  //             solver->apply_constraints();
  //             evaluation_point = local_evaluation_point;

  //             if (this->params.verbosity != Parameters::Verbosity::quiet)
  //               {
  //                 solver->pcout
  //                   << "\t\talpha value was kept at alpha = " << alpha
  //                   << " since alpha = " << alpha / 2
  //                   << " increased the residual" << std::endl;
  //               }
  //             current_res = last_alpha_res;
  //             break;
  //           }
  //         if (current_res < this->params.step_tolerance * last_res ||
  //             last_res < this->params.tolerance)
  //           {
  //             break;
  //           }
  //         last_alpha_res = current_res;
  //       }

  //     present_solution = evaluation_point;
  //     last_res         = current_res;
  //     ++outer_iteration;
  //   }
}

#endif
