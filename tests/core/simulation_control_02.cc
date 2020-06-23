/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

*
* Author: Bruno Blais, Polytechnique Montreal, 2020-
*/

// This test checks that the steady-state simulation control stops
// at the correct moment and behaves in a correct manner

#include <core/parameters.h>

#include "../tests.h"
#include "core/simulation_control.h"
#include "solvers/navier_stokes_solver_parameters.h"

int
main()
{
  try
    {
      initlog();

      Parameters::SimulationControl simulationControlParameters;

      simulationControlParameters.dt     = 0.01;
      simulationControlParameters.adapt  = false;
      simulationControlParameters.maxCFL = 99;
      simulationControlParameters.method =

        Parameters::SimulationControl::TimeSteppingMethod::bdf1;

      simulationControlParameters.timeEnd                = 999;
      simulationControlParameters.number_mesh_adaptation = 9;
      simulationControlParameters.output_name            = "test";
      simulationControlParameters.subdivision            = 7;
      simulationControlParameters.output_folder          = "canard";
      simulationControlParameters.output_frequency       = 8;

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
