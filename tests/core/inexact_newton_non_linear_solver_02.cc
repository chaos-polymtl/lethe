/**
 * @brief Test the solution of a non-linear problem, then resets
 * the solution vector of the problem to solve it again
 * while reusing the previous matrix.
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
    Parameters::NonLinearSolver::SolverType::inexact_newton,
    Parameters::NonLinearSolver::KinsolStrategy::
      normal_newton, // kinsol strategy, not used in this case
    1e-8,            // tolerance
    10,              // maxIter
    4,               // display precision
    false,           // force rhs calculation
    0.1,             // matrix tolerance
    0.99,            // step_tolerance
    true             // reuse matrix accross problems
  };

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<NonLinearProblemTestClass> solver =
    std::make_unique<NonLinearProblemTestClass>(params);


  deallog << "Solving non-linear system " << std::endl;
  // Solve the non-linear system of equation a first time
  solver->solve_non_linear_system(true);

  solver->reset();

  // Solve it again reusing the system matrix
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
