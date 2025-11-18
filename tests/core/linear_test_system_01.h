// SPDX-FileCopyrightText: Copyright (c) 2019-2021, 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II includes
#include <deal.II/lac/lapack_full_matrix.h>

// Lethe
#include <core/physics_solver.h>

/**
 * @brief The LinearProblemTestClass houses a very simple linear system that is used to test the disabled non-linear solver
 * The linear solution is obtained using LAPACK solve direct solver
 *
 * We use the base deal.II types (Vector<double>) to store the information as
 * there is no need to use the complex Trilinos Vector that support parallelism
 * in this case
 *
 * The system that is solved is :
 * x_0 + x_1 = 0
 *       2*x_1 = -3
 *
 */

class LinearProblemTestClass : public PhysicsSolver<Vector<double>>
{
public:
  LinearProblemTestClass(Parameters::NonLinearSolver &params)
    : PhysicsSolver(params)
  {
    // Initialize the vectors needed for the Physics Solver
    system_rhs.reinit(2);
    present_solution.reinit(2);
  }


  // Assembling Jacobian Matrix
  virtual void
  assemble_system_matrix() override
  {
    system_matrix.reinit(2);
    // System
    // x_0 +x_1 = 0
    // 2*x_1 = -3 
    //

    system_matrix.set(0, 0, 1);
    system_matrix.set(0, 1, 1);
    system_matrix.set(1, 0, 0);
    system_matrix.set(1, 1, 2);

    system_matrix.set_property(LAPACKSupport::general);
    system_matrix.compute_lu_factorization();
  }

  virtual void
  setup_preconditioner() override
  {}

  virtual void
  assemble_system_rhs() override
  {
    system_rhs[0] = 0;
    system_rhs[1] = -3;
  }

  /**
   * @brief solve_linear_system
   *
   * Solve the linear system of equation using LAPACK
   */

  void
  solve_linear_system(const bool) override
  {
    system_matrix.solve(system_rhs);
    present_solution = system_rhs;
  }

  virtual void
  apply_constraints() override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no constraints to apply.");
  }

  virtual Vector<double> &
  get_evaluation_point() override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no evaluation point available.");
  };
  virtual Vector<double> &
  get_local_evaluation_point() override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no local evaluation point available.");
  };
  virtual Vector<double> &
  get_newton_update() override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no newton update available.");
  };
  virtual void
  output_newton_update_norms(const unsigned int display_precision) override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no newton update norms available.");
  };
  virtual Vector<double> &
  get_present_solution() override
  {
    return present_solution;
  };
  virtual Vector<double> &
  get_system_rhs() override
  {
    return system_rhs;
  };
  virtual AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no constraints available.");
  };
  virtual double
  get_residual_rescale_metric() const override
  {
    throw std::runtime_error(
      "Newton iteration disabled, no residual rescale metric available.");
  }


private:
  LAPACKFullMatrix<double> system_matrix;
  Vector<double>           system_rhs;
  Vector<double>           present_solution;
};
