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
  Parameters::NonLinearSolver params{
    .verbosity = Parameters::Verbosity::verbose,
    .solver    = Parameters::NonLinearSolver::SolverType::disabled,
    .kinsol_strategy =
      Parameters::NonLinearSolver::KinsolStrategy::normal_newton, // not used in
                                                                  // this case
    .tolerance                    = 1e-8,   // not used in this case
    .max_iterations               = 10,     // not used in this case
    .display_precision            = 4,      // not used in this case
    .force_rhs_calculation        = false,  // not used in this case
    .matrix_tolerance             = 0.1,    // not used in this case
    .step_tolerance               = 0.99,   // not used in this case
    .reuse_matrix                 = false,  // not used in this case
    .reuse_preconditioner         = false,  // not used in this case
    .abort_at_convergence_failure = false}; // not used in this case

  deallog << "Creating solver" << std::endl;

  // Create an instantiation of the Test Class
  std::unique_ptr<LinearProblemTestClass> solver =
    std::make_unique<LinearProblemTestClass>(params);


  deallog << "Solving linear system " << std::endl;
  // We use solve_non_linear_system of equation even if the system is linear
  // because it inherits from the NonLinearSolver base class.
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
