// SPDX-FileCopyrightText: Copyright (c) 2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks that the time-step round-off errors
 * do not affect the end of the simulation
 */

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  Parameters::SimulationControl simulation_control_parameters;

  simulation_control_parameters.dt     = 0.02;
  simulation_control_parameters.adapt  = false;
  simulation_control_parameters.maxCFL = 1;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.time_step_independent_of_end_time = true;

  simulation_control_parameters.time_end                   = 20;
  simulation_control_parameters.number_mesh_adaptation     = 0;
  simulation_control_parameters.output_name                = "test";
  simulation_control_parameters.subdivision                = 7;
  simulation_control_parameters.output_folder              = "canard";
  simulation_control_parameters.output_iteration_frequency = 8;

  {
    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Constant time-stepping
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;
      }
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
