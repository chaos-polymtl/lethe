// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the capabilities of the SubSimulationControlDEM
 * to adaptively calculate the step in both the fixed coupling frequency
 * and in the fixed fraction of Rayleigh time step modes.
 */


// Lethe
#include <core/sub_simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  // We first create a case with fixed number of iterations
  SubSimulationControlDEM fixed_iterations(
    SubSimulationControlDEM::DEMSubIterationLogic::fixed_number_of_iterations,
    0.1,
    2,
    1,
    1);

  deallog << "Fixed iteration number" << std::endl;

  while (fixed_iterations.iterate())
    {
      deallog << "Iteration: " << fixed_iterations.get_iteration()
              << " Time step: " << fixed_iterations.get_time_step()
              << std::endl;
      ;
    }

  // We then create a case with a fixed fraction of the Rayleigh time step
  SubSimulationControlDEM fixed_fraction(
    SubSimulationControlDEM::DEMSubIterationLogic::
      fixed_fraction_of_rayleigh_time_step,
    0.1,
    11,
    0.1,
    0.6);

  deallog << "Fixed fraction of Rayleigh time step" << std::endl;

  while (fixed_fraction.iterate())
    {
      deallog << "Iteration: " << fixed_fraction.get_iteration()
              << " Time step: " << fixed_fraction.get_time_step() << std::endl;
    }
}


int
main()
{
  try
    {
      initlog();
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
}
