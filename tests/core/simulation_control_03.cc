// SPDX-FileCopyrightText: Copyright (c) 2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks whether the is_output_iteration behaves correctly for iteration output control and time output control (in the case of constant time step and adaptive time step).
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

  simulation_control_parameters.dt     = 0.1;
  simulation_control_parameters.adapt  = false;
  simulation_control_parameters.maxCFL = 2;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.time_end = 0.5;
  simulation_control_parameters.output_control =
    Parameters::SimulationControl::OutputControl::iteration;
  simulation_control_parameters.output_iteration_frequency        = 2;
  simulation_control_parameters.time_step_independent_of_end_time = true;
  simulation_control_parameters.output_time_interval              = {0,
                                                                     1.7976931348623157e3};

  {
    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - constant output
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
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {0.3};
    simulation_control_parameters.output_time_frequency = -1;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - specific time output
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - specific time output" << std::endl;
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
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = -1;
    simulation_control_parameters.output_time_interval  = {0.2, 0.4};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - time interval output
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - specific time interval" << std::endl;
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
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = 0.2;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - time interval output
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - output time frequency" << std::endl;
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
    simulation_control_parameters.time_end                     = 20;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.max_dt                       = 1e6;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::iteration;
    simulation_control_parameters.output_iteration_frequency = 8;
    simulation_control_parameters.output_time_interval       = {0,
                                                                1.7976931348623157e3};

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
    simulation_control_parameters.time_end                     = 25;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {7.5, 22.1};
    simulation_control_parameters.output_time_frequency = -1;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Adaptative time-stepping - specific time output
    deallog << "*************************************************" << std::endl;
    deallog << "Adaptative time stepping - specific time output" << std::endl;
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
    simulation_control_parameters.time_end                     = 25;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.max_dt                       = 1e6;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = -1;
    simulation_control_parameters.output_time_interval  = {7.5, 17};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Adaptative time-stepping - time interval output
    deallog << "*************************************************" << std::endl;
    deallog << "Adaptative time stepping - interval time output" << std::endl;
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
    simulation_control_parameters.time_end                     = 25;
    simulation_control_parameters.dt                           = 1;
    simulation_control_parameters.adaptative_time_step_scaling = 1.2;
    simulation_control_parameters.maxCFL                       = 2;
    simulation_control_parameters.max_dt                       = 1e6;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = 4;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Adaptative time-stepping - output time frequency
    deallog << "*************************************************" << std::endl;
    deallog << "Adaptative time stepping - output time frequency" << std::endl;
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
