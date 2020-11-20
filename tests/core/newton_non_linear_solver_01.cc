#include "non_linear_test_system_01.h"

/**
 * @brief The TestClass tests the non-linear solvers using a simple systme of two
 * equations, only one of which is non-linear
 */

void
test(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  Parameters::NonLinearSolver params{
    Parameters::Verbosity::quiet,
    Parameters::NonLinearSolver::SolverType::newton,
    1e-8, // tolerance
    10,   // maxIter
    4,    // display precision
    1     // skip iterations
  };

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<TestClass> solver = std::make_unique<TestClass>(params);


  deallog << "Solving non-linear system " << std::endl;
  // Solve the non-linear system of equation
  solver->solve_non_linear_system(
    Parameters::SimulationControl::TimeSteppingMethod::steady, true, true);

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
      test(argc, argv);
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
