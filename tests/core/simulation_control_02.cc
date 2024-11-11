// SPDX-FileCopyrightText: Copyright (c) 2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks that the steady-state simulation control stops
 * at the correct moment and behaves in a correct manner
 */

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  Parameters::SimulationControl simulationControlParameters;

  simulationControlParameters.dt     = 0.01;
  simulationControlParameters.adapt  = false;
  simulationControlParameters.maxCFL = 99;
  simulationControlParameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;

  simulationControlParameters.time_end                   = 999;
  simulationControlParameters.number_mesh_adaptation     = 9;
  simulationControlParameters.output_name                = "test";
  simulationControlParameters.subdivision                = 7;
  simulationControlParameters.output_folder              = "canard";
  simulationControlParameters.output_iteration_frequency = 8;
  simulationControlParameters.output_time_interval       = {0, 1000000000};
  simulationControlParameters.time_step_independent_of_end_time = true;

  SimulationControlSteady simulation_control(simulationControlParameters);

  deallog << "Iteration : " << simulation_control.get_step_number()
          << std::endl;

  while (simulation_control.integrate())
    {
      deallog << "Iteration : " << simulation_control.get_step_number()
              << "    Time : " << simulation_control.get_current_time()
              << std::endl;
      if (simulation_control.is_at_start())
        deallog << "This is the first iteration" << std::endl;

      if (simulation_control.is_output_iteration())
        deallog << "This is an output iteration" << std::endl;
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
