// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks that a transient simulation ends at the right moment
 * for both end control types. With the iteration end control, the simulation
 * must end after the prescribed number of transient iterations and the end time
 * must have no effect whatsoever, including on the value of the time step. With
 * the time end control, the simulation must end at the end time, which is also
 * the behavior obtained when no end control is specified.
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

  simulation_control_parameters.dt                            = 0.1;
  simulation_control_parameters.time_step_adaptation_required = false;
  simulation_control_parameters.adapt_with_cfl                = false;
  simulation_control_parameters.maxCFL                        = 2;
  simulation_control_parameters.max_dt                        = 1e6;
  simulation_control_parameters.adaptative_time_step_scaling  = 1.2;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.startup_timestep_scaling = 0.4;
  simulation_control_parameters.bdf_startup_method =
    Parameters::SimulationControl::BDFStartupMethods::multiple_step_bdf;
  simulation_control_parameters.time_end                             = 0.5;
  simulation_control_parameters.number_mesh_adaptation               = 0;
  simulation_control_parameters.time_step_independent_of_end_time    = true;
  simulation_control_parameters.adapt_with_capillary_time_step_ratio = false;

  {
    // The end control is not specified. The simulation must end at the end
    // time, as it did before the end control parameter was introduced.
    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Default end control - time end at 0.5" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;
      }
  }

  {
    // The simulation ends at a given iteration. The end time is larger than the
    // time reached at the last iteration and must not end the simulation.
    simulation_control_parameters.end_control =
      Parameters::SimulationControl::EndControl::iteration;
    simulation_control_parameters.iteration_end = 3;
    simulation_control_parameters.time_end      = 100;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Iteration end control - iteration end at 3" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;
      }
  }

  {
    // The simulation ends at a given iteration, but the end time is now reached
    // before the last iteration. The end time must still have no effect.
    simulation_control_parameters.iteration_end = 8;
    simulation_control_parameters.time_end      = 0.25;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Iteration end control - time end reached before the last "
               "iteration"
            << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;
      }
  }

  {
    // Same case, but the time step is not independent of the end time. Since
    // the simulation ends at a given iteration, the time step must not be
    // clipped to reach the end time exactly.
    simulation_control_parameters.time_step_independent_of_end_time = false;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Iteration end control - time step not independent of end time"
            << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;
      }
  }

  {
    // Adaptive time-stepping with an iteration end control. The number of
    // iterations is unaffected by the growth of the time step.
    simulation_control_parameters.time_step_independent_of_end_time = true;
    simulation_control_parameters.time_step_adaptation_required     = true;
    simulation_control_parameters.adapt_with_cfl                    = true;
    simulation_control_parameters.iteration_end                     = 5;
    simulation_control_parameters.dt                                = 1;
    simulation_control_parameters.time_end                          = 3;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Iteration end control - adaptative time stepping" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << "    Time step : " << simulation_control.get_time_step()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << "    Time step : " << simulation_control.get_time_step()
                << std::endl;

        simulation_control.set_CFL(simulation_control.get_time_step());
      }
  }

  {
    // An iteration end control set to zero iteration ends the simulation
    // before it starts.
    simulation_control_parameters.time_step_adaptation_required = false;
    simulation_control_parameters.adapt_with_cfl                = false;
    simulation_control_parameters.iteration_end                 = 0;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Iteration end control - iteration end at 0" << std::endl;
    deallog << "*************************************************" << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
                << "    Time : " << simulation_control.get_current_time()
                << std::endl;
      }

    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;
  }

  {
    // Switching back to the time end control ends the simulation at the end
    // time, no matter the value of the end iteration.
    simulation_control_parameters.end_control =
      Parameters::SimulationControl::EndControl::time;
    simulation_control_parameters.iteration_end = 2;
    simulation_control_parameters.dt            = 0.1;
    simulation_control_parameters.time_end      = 0.5;

    SimulationControlTransient simulation_control(
      simulation_control_parameters);

    deallog << "*************************************************" << std::endl;
    deallog << "Time end control - iteration end has no effect" << std::endl;
    deallog << "*************************************************" << std::endl;
    deallog << "Iteration : " << simulation_control.get_iteration_number()
            << "    Time : " << simulation_control.get_current_time()
            << std::endl;

    while (simulation_control.integrate())
      {
        deallog << "Iteration : " << simulation_control.get_iteration_number()
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
