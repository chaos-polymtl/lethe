// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief The TestClass tests the non-linear solvers using a simple system of two
 * equations, only one of which is non-linear
 */

// Lethe
#include <core/parameters.h>

// Tests
#include <../tests/core/linear_test_system_01.h>

#include <../tests/tests.h>

void
test()
{

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<LinearProblemTestClass> solver =
    std::make_unique<LinearProblemTestClass>();


  deallog << "Solving linear system " << std::endl;
  solver->solve_governing_system();

  auto &solution = solver->get_system_rhs();
  deallog << "The final solution is : " << solution[0] << " "
          << solution[1] << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
