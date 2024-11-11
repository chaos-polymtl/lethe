// SPDX-FileCopyrightText: Copyright (c) 2019-2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the read and write capacities of the
 * SimulationFlowControl class
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
  simulationControlParameters.time_end                          = 999;
  simulationControlParameters.number_mesh_adaptation            = 9;
  simulationControlParameters.output_name                       = "test";
  simulationControlParameters.subdivision                       = 7;
  simulationControlParameters.output_folder                     = "canard";
  simulationControlParameters.output_iteration_frequency        = 8;
  simulationControlParameters.time_step_independent_of_end_time = true;

  SimulationControlTransient simulationControl(simulationControlParameters);

  for (int i = 0; i < 10; ++i)
    simulationControl.integrate();

  simulationControl.save("testFile");
  simulationControl.read("testFile");

  deallog << "dt                  : " << simulationControl.get_time_step()
          << std::endl;
  deallog << "CFL                 : " << simulationControl.get_CFL()
          << std::endl;
  deallog << "time                : " << simulationControl.get_current_time()
          << std::endl;
  deallog << "iter                : " << simulationControl.get_step_number()
          << std::endl;
  deallog << "OK" << std::endl;
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
