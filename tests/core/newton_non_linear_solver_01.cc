#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/vector.h>

#include <core/newton_non_linear_solver.h>
#include <core/physics_solver.h>
#include <core/simulation_control.h>

#include <iostream>
#include <memory>

#include "../tests.h"

/**
 * @brief The TestClass tests the non-linear solvers using a simple systme of two
 * equations, only one of which is non-linear
 * It uses a FullMatrix and a Vector to store the analytical jacobian et the
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
    : PhysicsSolver(new NewtonNonLinearSolver<Vector<double>>(this, params))
  {
    // Initialize the vectors needed for the Physics Solver
    this->evaluation_point.reinit(2);
    this->system_rhs.reinit(2);
    this->local_evaluation_point.reinit(2);
    this->present_solution.reinit(2);
    this->newton_update.reinit(2);

    // Set the initial value of the solution
    this->present_solution[0] = 1;
    this->present_solution[1] = 0;

    this->nonzero_constraints.close();
  }


  // Assembling Jacobian Matrix
  virtual void
  assemble_matrix_and_rhs(
    const Parameters::SimulationControl::TimeSteppingMethod
      time_stepping_method) override
  {
    system_matrix.reinit(2);
    // System
    // x_0*x_0 +x_1 = 0
    // 2*x_1 + 3 = 0
    //
    // Jacobian
    // 2x_0     1
    // 0        2
    double x_0 = this->evaluation_point[0];
    double x_1 = this->evaluation_point[1];
    system_matrix.set(0, 0, 2 * x_0);
    system_matrix.set(0, 1, 1);
    system_matrix.set(1, 0, 0);
    system_matrix.set(1, 1, 2);

    this->system_rhs[0] = -(x_0 * x_0 + x_1);
    this->system_rhs[1] = -(2 * x_1 + 3);
  }

  virtual void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override
  {
    double x_0          = this->evaluation_point[0];
    double x_1          = this->evaluation_point[1];
    this->system_rhs[0] = -(x_0 * x_0 + x_1);
    this->system_rhs[1] = -(2 * x_1 + 3);
  }

  /**
   * @brief solve_linear_system
   *
   * Solve the linear system of equation using LAPACK
   */

  void
  solve_linear_system(const bool, const bool) override
  {
    system_matrix.set_property(LAPACKSupport::general);
    //    system_matrix.compute_lu_factorization();
    system_matrix.compute_lu_factorization();
    //    system_matrix.apply_lu_factorization(this->system);
    system_matrix.solve(this->system_rhs);
    newton_update = this->system_rhs;
  }

  virtual void
  apply_constraints()
  {}

private:
  LAPACKFullMatrix<double> system_matrix;
  Vector<double>           rhs;
};

int
main()
{
  initlog();

  Parameters::NonLinearSolver params{
    Parameters::Verbosity::quiet,
    Parameters::NonLinearSolver::SolverType::newton,
    1e-6, // tolerance
    10,   // maxIter
    4     // display precision
  };

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<TestClass> solver = std::make_unique<TestClass>(params);


  deallog << "Solving non-linear system " << std::endl;
  // Solve the non-linear system of equation
  solver->solve_non_linear_system(
    Parameters::SimulationControl::TimeSteppingMethod::steady, true);


  deallog << "The final solution is : " << solver->get_present_solution()
          << std::endl;
}
