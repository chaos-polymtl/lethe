#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/vector.h>

#include <core/basic_non_linear_solver.h>
#include <core/physics_solver.h>
#include <core/simulationcontrol.h>

#include <iostream>
#include <memory>

#include "../tests.h"

class TestClass : public PhysicsSolver<TrilinosWrappers::MPI::Vector>
{
public:
  TestClass(Parameters::NonLinearSolver &params)
    : PhysicsSolver(new BasicNonLinearSolver(this, params, 1e-15, 1e-15))
  {}

  void
  assemble_matrix_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                        time_stepping_method) override
  {
    // Assemblage de la jacobienne

    // Systeme:
    // x_0*x_0 +x_1 = 0
    // 2*x_1 + 3 = 0
    // J = [2x0, 1 \n 0, 2]
    // x_init = [1 \n 3]
    // J(x_init) = [2, 1 \n 0, 2]

    // Pseudocode-ish:
    // system_rhs = FullMatrix(2, 2);
    // system_rhs[0][0] = [2];
    // system_rhs[0][1] = [1];
    // system_rhs[1][0] = [0];
    // system_rhs[1][1] = [2];
  }

  void
  assemble_rhs(const Parameters::SimulationControl::TimeSteppingMethod
                 time_stepping_method) override
  {
    // Est-ce que ici on ferait evaluation_point = [1, 3]?
  }

  void
  solve_linear_system(const bool   initial_step,
                      const double absolute_residual,
                      const double relative_residual) override
  {
    // Est-ce que ici on utiliserait present_solution?
  }
};

int
main()
{
  std::cout << "bonjour" << std::endl;
  Parameters::NonLinearSolver params{
    Parameters::Verbosity::quiet,
    0.1, // tolerance
    10,  // maxIter
    4    // display precision
  };
  std::unique_ptr<PhysicsSolver<TrilinosWrappers::MPI::Vector>> solver =
    std::make_unique<TestClass>(TestClass(params));
  solver->solve_non_linear_system(
    Parameters::SimulationControl::TimeSteppingMethod::steady, true);

  std::cout << "fini!" << std::endl;
}