// SPDX-FileCopyrightText: Copyright (c) 2020, 2024, 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks whether the is_output_iteration behaves correctly for iteration output control and time output control (in the case of constant time step and adaptive time step). The output is also checked at the initial iteration. For the output time frequency, the test also checks that the initial and last iterations are written, that the output times do not drift when the output time frequency is not a multiple of the time step, that the output time interval is respected, and that repeated calls within the same iteration return the same result. For specific output times, the test also checks that unsorted output times, an output time at the initial time, and several output times within the same time step are handled, and that repeated calls within the same iteration return the same result.
 */

// Lethe
#include <core/parameters.h>
#include <core/simulation_control.h>

// deal.II
#include <deal.II/base/exceptions.h>

// Tests (with common definitions)
#include <../tests/tests.h>

/**
 * @brief Query whether the current iteration is an output iteration twice and
 * check that both queries return the same result.
 *
 * @param[in,out] simulation_control Simulation control that is queried.
 *
 * @return true if the current iteration is an output iteration.
 */
bool
is_output_iteration_queried_twice(
  SimulationControlTransient &simulation_control)
{
  const bool first_query  = simulation_control.is_output_iteration();
  const bool second_query = simulation_control.is_output_iteration();
  AssertThrow(first_query == second_query,
              ExcMessage("Two calls to is_output_iteration within the same "
                         "iteration returned different results."));
  return first_query;
}

void
test()
{
  Parameters::SimulationControl simulation_control_parameters;

  simulation_control_parameters.dt                            = 0.1;
  simulation_control_parameters.time_step_adaptation_required = false;
  simulation_control_parameters.adapt_with_cfl                = false;
  simulation_control_parameters.maxCFL                        = 2;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.end_control =
    Parameters::SimulationControl::EndControl::time;
  simulation_control_parameters.time_end = 0.5;
  simulation_control_parameters.output_control =
    Parameters::SimulationControl::OutputControl::iteration;
  simulation_control_parameters.output_iteration_frequency           = 2;
  simulation_control_parameters.time_step_independent_of_end_time    = true;
  simulation_control_parameters.output_time_interval                 = {0,
                                                                        1.7976931348623157e3};
  simulation_control_parameters.adapt_with_capillary_time_step_ratio = false;

  {
    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - constant output
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - constant output" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (simulation_control.is_output_iteration())
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (simulation_control.is_output_iteration())
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
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

    // Constant time-stepping - output time frequency. The initial and last
    // iterations must be output iterations.
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - output time frequency" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;
      }
  }

  {
    simulation_control_parameters.time_step_adaptation_required = true;
    simulation_control_parameters.adapt_with_cfl                = true;
    simulation_control_parameters.time_end                      = 20;
    simulation_control_parameters.dt                            = 1;
    simulation_control_parameters.adaptative_time_step_scaling  = 1.2;
    simulation_control_parameters.maxCFL                        = 2;
    simulation_control_parameters.max_dt                        = 1e6;
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    if (simulation_control.is_output_iteration())
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
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
    simulation_control_parameters.time_step_adaptation_required = true;
    simulation_control_parameters.adapt_with_cfl                = true;
    simulation_control_parameters.time_end                      = 25;
    simulation_control_parameters.dt                            = 1;
    simulation_control_parameters.adaptative_time_step_scaling  = 1.2;
    simulation_control_parameters.maxCFL                        = 2;
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;

        simulation_control.set_CFL(simulation_control.get_time_step());
      }
  }

  {
    simulation_control_parameters.time_step_adaptation_required = true;
    simulation_control_parameters.adapt_with_cfl                = true;
    simulation_control_parameters.time_end                      = 25;
    simulation_control_parameters.dt                            = 1;
    simulation_control_parameters.adaptative_time_step_scaling  = 1.2;
    simulation_control_parameters.maxCFL                        = 2;
    simulation_control_parameters.max_dt                        = 1e6;
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    if (simulation_control.is_output_iteration())
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
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
    simulation_control_parameters.time_step_adaptation_required = true;
    simulation_control_parameters.adapt_with_cfl                = true;
    simulation_control_parameters.time_end                      = 25;
    simulation_control_parameters.dt                            = 1;
    simulation_control_parameters.adaptative_time_step_scaling  = 1.2;
    simulation_control_parameters.maxCFL                        = 2;
    simulation_control_parameters.max_dt                        = 1e6;
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
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;

        simulation_control.set_CFL(simulation_control.get_time_step());
      }
  }

  {
    simulation_control_parameters.time_step_adaptation_required     = false;
    simulation_control_parameters.adapt_with_cfl                    = false;
    simulation_control_parameters.time_step_independent_of_end_time = false;
    simulation_control_parameters.time_end                          = 1.0;
    simulation_control_parameters.dt                                = 0.1;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = 0.25;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - output time frequency that is not a multiple of
    // the time step. The output times (0.25, 0.5, 0.75 and 1) must not drift
    // and the end time, which is also an output time, must be written.
    deallog << "*************************************************" << std::endl;
    deallog
      << "Constant time stepping - output time frequency not multiple of time step"
      << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;
      }
  }

  {
    simulation_control_parameters.time_step_adaptation_required     = false;
    simulation_control_parameters.adapt_with_cfl                    = false;
    simulation_control_parameters.time_step_independent_of_end_time = false;
    simulation_control_parameters.time_end                          = 1.0;
    simulation_control_parameters.dt                                = 0.1;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {-1};
    simulation_control_parameters.output_time_frequency = 0.2;
    simulation_control_parameters.output_time_interval  = {0.5, 0.8};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - output time frequency within an output time
    // interval. The output times start at the lower bound of the interval
    // (0.5 and 0.7), and the initial and last iterations are not written since
    // they lie outside of the interval.
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - output time frequency within interval"
            << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;
      }
  }

  {
    simulation_control_parameters.time_step_adaptation_required     = false;
    simulation_control_parameters.adapt_with_cfl                    = false;
    simulation_control_parameters.time_step_independent_of_end_time = false;
    simulation_control_parameters.time_end                          = 1.0;
    simulation_control_parameters.dt                                = 0.1;
    simulation_control_parameters.output_control =
      Parameters::SimulationControl::OutputControl::time;
    simulation_control_parameters.output_times_vector   = {0.75, 0, 0.32, 0.31};
    simulation_control_parameters.output_time_frequency = -1;
    simulation_control_parameters.output_time_interval  = {0,
                                                           1.7976931348623157e3};

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    // Constant time-stepping - specific output times that are not sorted,
    // include the initial time, and of which two (0.31 and 0.32) lie within the
    // same time step. The initial iteration must be written, as well as the
    // iterations before and after each other output time.
    deallog << "*************************************************" << std::endl;
    deallog << "Constant time stepping - specific output times" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    if (is_output_iteration_queried_twice(simulation_control))
      deallog << "This is an output iteration" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;

        if (simulation_control.is_at_start())
          deallog << "This is the first time step" << std::endl;

        if (is_output_iteration_queried_twice(simulation_control))
          deallog << "This is an output iteration" << std::endl;
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
