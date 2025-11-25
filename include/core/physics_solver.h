// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_physics_solver_h
#define lethe_physics_solver_h


#include <core/inexact_newton_non_linear_solver_strategy.h>
#include <core/kinsol_newton_non_linear_solver_strategy.h>
#include <core/linear_solver_strategy.h>
#include <core/newton_non_linear_solver_strategy.h>
#include <core/parameters.h>
#include <core/physics_solver_strategy.h>

#include <deal.II/lac/affine_constraints.h>

/**
 * @brief Contain all the common elements of a physics solver. It creates
 * the non-linear solver as specified by the user using the parameters file and
 * provides all the necessary elements needed by the solver to solve a physics
 * problem.
 *
 * @param[in] non_linear_solver_parameters Set of parameters that will be used
 * to construct the non-linear solver.
 */
template <typename VectorType>
class PhysicsSolver
{
public:
  /**
   * @brief Constructor for the non-linear physics.
   *
   * @param[in] Non-linear solver parameters as specified in the
   * simulation parameter file.
   *
   */
  PhysicsSolver(const Parameters::NonLinearSolver non_linear_solver_parameters);


  /**
   * @brief Constructor for the linear physics. Since the Physics is linear, a LinearSolutionStrategy is automatically generated as a PhysicsSolverStrategy.
   *
   */
  PhysicsSolver();


  /**
   * @brief Destructor.
   *
   */
  virtual ~PhysicsSolver()
  {
    delete physics_solving_strategy;
  }

  /**
   * @brief Assemble the matrix.
   */
  virtual void
  assemble_system_matrix() = 0;

  /**
   * @brief Assemble the right-hand side (rhs).
   */
  virtual void
  assemble_system_rhs() = 0;

  /**
   * @brief Set up the preconditioner.
   *
   */
  virtual void
  setup_preconditioner() = 0;

  /**
   * @brief Solve the linear system of equations.
   */
  virtual void
  solve_linear_system() = 0;

  /**
   * @brief Solve the global system of equations according to a given strategy, either linear or not.
   */
  void
  solve_governing_system();

  /**
   * @brief Applies constraints to a local_evaluation_point.
   */
  virtual void
  apply_constraints()
  {
    auto &nonzero_constraints    = get_nonzero_constraints();
    auto &local_evaluation_point = get_local_evaluation_point();
    nonzero_constraints.distribute(local_evaluation_point);
  }

  /**
   * @brief Getter methods that give access to the private attributes of the
   * physics being solved. These methods must be provided by the physics as they
   * are the key to solve a problem using Lethe.
   */
  virtual VectorType &
  get_evaluation_point() = 0;
  virtual VectorType &
  get_local_evaluation_point() = 0;
  virtual VectorType &
  get_newton_update() = 0;
  virtual VectorType &
  get_present_solution() = 0;
  virtual VectorType &
  get_system_rhs() = 0;
  virtual AffineConstraints<double> &
  get_nonzero_constraints() = 0;

  /**
   * @brief Output the L2 and Linfty norms of the correction vector.
   *
   * @param[in] display_precision Number of outputted digits.
   */
  virtual void
  output_newton_update_norms(const unsigned int display_precision) = 0;

  /**
   * @brief Return the metric for residual rescaling. By default, should return 1.
   * If the rescale_residual_by_volume is set to true, the method
   * returns the square root of the global volume of the triangulation.
   *
   * @return Rescale metric.
   */
  virtual double
  get_residual_rescale_metric() const = 0;

  /**
   * @brief Default way to evaluate the residual for the nonlinear solver, which is simply the norm of the RHS vector.
   * Some application may use more complex evaluation of the residual and
   * override this method.
   *
   * @return l2 norm of the RHS vector.
   */
  virtual double
  get_current_residual()
  {
    auto &system_rhs = get_system_rhs();
    return system_rhs.l2_norm();
  }

  /**
   * @brief Return the current newton iteration of this physics solver.
   */
  inline unsigned int
  get_current_newton_iteration() const
  {
    return physics_solving_strategy->get_current_newton_iteration();
  }

  ConditionalOStream                                pcout;
  Parameters::SimulationControl::TimeSteppingMethod time_stepping_method;

private:
  PhysicsSolverStrategy<VectorType> *physics_solving_strategy;
};

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver(
  const Parameters::NonLinearSolver non_linear_solver_parameters)
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
{
  switch (non_linear_solver_parameters.solver)
    {
      case Parameters::NonLinearSolver::SolverType::newton:
        physics_solving_strategy =
          new NewtonNonLinearSolverStrategy<VectorType>(
            this, non_linear_solver_parameters);
        break;
      case Parameters::NonLinearSolver::SolverType::kinsol_newton:
        physics_solving_strategy =
          new KinsolNewtonNonLinearSolverStrategy<VectorType>(
            this, non_linear_solver_parameters);
        break;
      case Parameters::NonLinearSolver::SolverType::inexact_newton:
        physics_solving_strategy =
          new InexactNewtonNonLinearSolverStrategy<VectorType>(
            this, non_linear_solver_parameters);
        break;
      default:
        break;
    }
}

template <typename VectorType>
PhysicsSolver<VectorType>::PhysicsSolver()
  : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
{
  physics_solving_strategy = new LinearSolverStrategy<VectorType>(this);
}


template <typename VectorType>
void
PhysicsSolver<VectorType>::solve_governing_system()
{
  {
    this->physics_solving_strategy->solve();
  }
}
#endif
