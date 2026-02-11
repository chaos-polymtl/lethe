// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the capabilities of the SubSimulationControlDEM
 * to correctly calculate the time step and the number of iterations in both
 * the fixed coupling frequency and in the fixed fraction of Rayleigh time step
 * modes. In both cases, the test iterates through all DEM sub iterations and
 * verifies the iteration count and time step value at each step.
 */


// Lethe
#include <core/sub_simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  // Test 1: Fixed number of iterations mode
  // With a time interval of 0.1 s and a coupling frequency of 2, the DEM
  // time step should be 0.1 / 2 = 0.05 s and 2 iterations should be performed.
  // The Rayleigh time (1 s) and its fraction (1) are unused in this mode.
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

  // Test 2: Fixed fraction of Rayleigh time step mode
  // With a Rayleigh characteristic time of 0.1 s and a target fraction of 0.6,
  // the DEM time step should be 0.1 * 0.6 = 0.06 s. With a time interval of
  // 0.1 s, the number of iterations is ceil(0.1 / 0.06) = 2, and the actual
  // time step is adjusted to 0.1 / 2 = 0.05 s to evenly divide the interval.
  // The coupling frequency (11) is unused in this mode.
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
