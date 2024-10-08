// SPDX-FileCopyrightText: Copyright (c) 2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks that the transient simulation control stops
 * at the correct moment and behaves in a correct manner.
 */

// We first check if the constant time-stepping number work
// then we test adaptative time-stepping

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>

void
test()
{
  Parameters::SimulationControl simulation_control_parameters;

  simulation_control_parameters.dt     = 0.1;
  simulation_control_parameters.adapt  = false;
  simulation_control_parameters.maxCFL = 2;
  simulation_control_parameters.method =

    Parameters::SimulationControl::TimeSteppingMethod::bdf1;

  simulation_control_parameters.timeEnd                = 0.5;
  simulation_control_parameters.number_mesh_adaptation = 9;
  simulation_control_parameters.output_name            = "test";
  simulation_control_parameters.subdivision            = 7;
  simulation_control_parameters.output_folder          = "canard";
  simulation_control_parameters.output_frequency       = 8;
  simulation_control_parameters.output_time_interval   = {0, 1000000000};

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

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (simulation_control.is_output_iteration())
          deallog << "This is an output iteration" << std::endl;
      }
  }

  {
    simulation_control_parameters.adapt                        = true;
    simulation_control_parameters.timeEnd                      = 20;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.max_dt                       = 1e6;


    SimulationControlTransient simulation_control(
      simulation_control_parameters);


    // Adaptative time-stepping - constant output
    deallog << "*************************************************" << std::endl;
    deallog << "Adaptative time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (simulation_control.is_output_iteration())
          deallog << "This is an output iteration" << std::endl;

        simulation_control.set_CFL(simulation_control.get_time_step());
      }
  }

  {
    simulation_control_parameters.adapt                        = true;
    simulation_control_parameters.timeEnd                      = 30;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.output_time                  = 7.5;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.max_dt = 1e6;


    SimulationControlTransientDynamicOutput simulation_control(
      simulation_control_parameters);


    // Adaptative time-stepping - time output
    deallog << "*************************************************" << std::endl;
    deallog << "Adaptative time stepping - time output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_step_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_step_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (simulation_control.is_output_iteration())
          deallog << "This is an output iteration" << std::endl;

        simulation_control.set_CFL(simulation_control.get_time_step());
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
