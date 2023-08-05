/**
 * @brief The TestClass tests the non-linear solvers using a simple system of two
 * equations, only one of which is non-linear
 */

// Lethe
#include <core/parameters.h>

// Tests (with common definitions)
#include <../tests/core/non_linear_test_system_01.h>

#include <../tests/tests.h>

void
test()
{
  Parameters::NonLinearSolver params{
    Parameters::Verbosity::verbose,
    Parameters::NonLinearSolver::SolverType::newton,
    Parameters::NonLinearSolver::KinsolStrategy::
      normal_newton, // kinsol strategy, not used in this case
    1e-8,            // tolerance
    10,              // maxIter
    4,               // display precision
    false,           // force rhs calculation
    0.1,             // matrix tolerance
    0.99,            // step_tolerance
    false            // reuse matrix accross problems

  };

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<NonLinearProblemTestClass> solver =
    std::make_unique<NonLinearProblemTestClass>(params);


  deallog << "Solving non-linear system " << std::endl;
  // Solve the non-linear system of equation
  solver->solve_non_linear_system(true);

  auto &present_solution = solver->get_present_solution();
  deallog << "The final solution is : " << present_solution[0] << " "
          << present_solution[1] << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
