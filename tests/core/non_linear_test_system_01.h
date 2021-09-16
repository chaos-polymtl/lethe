/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

 *
 * Author: Simon Gauvin, Polytechnique Montreal, 2019
 */

// Deal.II includes
#include <deal.II/lac/lapack_full_matrix.h>

// Lethe
#include <core/physics_solver.h>

/**
 * @brief The TestClass houses a very simple non-linear system that is used to test the various non-linear solvers
 * The linear solution is obtained using LAPACK
 * it uses a LAPACKMatrix and a Vector to store the analytical jacobian et the
 * right-hand side respectively
 *
 * We use the base deal.II types (Vector<double>) to store the information as
 * there is no need to use the complex Trilinos Vector that support parallelism
 * in this case
 *
 * The system that is solved is :
 * x_0^2 + x_1 = 0
 *       2*x_1 + 3 = 0
 *
 */

class TestClass : public PhysicsSolver<Vector<double>>
{
public:
  TestClass(Parameters::NonLinearSolver &params)
    : PhysicsSolver(params)
  {
    // Initialize the vectors needed for the Physics Solver
    evaluation_point.reinit(2);
    system_rhs.reinit(2);
    local_evaluation_point.reinit(2);
    present_solution.reinit(2);
    newton_update.reinit(2);

    // Set the initial value of the solution
    present_solution[0] = 1;
    present_solution[1] = 0;
  }


  // Assembling Jacobian Matrix
  virtual void
  assemble_system_matrix() override
  {
    system_matrix.reinit(2);
    // System
    // x_0*x_0 +x_1 = 0
    // 2*x_1 + 3 = 0
    //
    // Jacobian
    // 2x_0     1
    // 0        2
    double x_0 = evaluation_point[0];
    double x_1 = evaluation_point[1];
    system_matrix.set(0, 0, 2 * x_0);
    system_matrix.set(0, 1, 1);
    system_matrix.set(1, 0, 0);
    system_matrix.set(1, 1, 2);

    system_rhs[0] = -(x_0 * x_0 + x_1);
    system_rhs[1] = -(2 * x_1 + 3);

    system_matrix.set_property(LAPACKSupport::general);
    system_matrix.compute_lu_factorization();
  }

  virtual void
  assemble_system_rhs() override
  {
    double x_0 = evaluation_point[0];
    double x_1 = evaluation_point[1];

    system_rhs[0] = -(x_0 * x_0 + x_1);
    system_rhs[1] = -(2 * x_1 + 3);
  }

  /**
   * @brief solve_linear_system
   *
   * Solve the linear system of equation using LAPACK
   */

  void
  solve_linear_system(const bool, const bool) override
  {
    system_matrix.solve(system_rhs);
    newton_update = system_rhs;
  }

  virtual void
  apply_constraints()
  {}

  virtual Vector<double> &
  get_evaluation_point() override
  {
    return evaluation_point;
  };
  virtual Vector<double> &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  };
  virtual Vector<double> &
  get_newton_update() override
  {
    return newton_update;
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
    return dummy_constraints;
  };


private:
  LAPACKFullMatrix<double>  system_matrix;
  AffineConstraints<double> dummy_constraints;
  Vector<double>            system_rhs;
  Vector<double>            newton_update;
  Vector<double>            evaluation_point;
  Vector<double>            local_evaluation_point;
  Vector<double>            present_solution;
};
