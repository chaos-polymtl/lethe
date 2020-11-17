/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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
* Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2020-
*/

/**
 * @brief This code tests the reynolds stress calculations in 3d with
 * Trilinos vectors.
 */

#include <core/parameters.h>
#include <solvers/postprocessing_velocities.h>

#include "../tests.h"

void
test(int argc, char **argv)
{
  // MPI initialization
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  // SimulationControl parameters
  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.dt               = 0.1;
  simulation_control_parameters.timeEnd          = 1.0;
  simulation_control_parameters.adapt            = false;
  simulation_control_parameters.output_frequency = 1;

  // Variables for AverageVelocities
  AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet>
    postprocessing_velocities;

  auto simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  IndexSet locally_owned_dofs;
  locally_owned_dofs.add_range(0, 8);

  Parameters::PostProcessing postprocessing_parameters;
  postprocessing_parameters.calculate_average_velocities = true;
  postprocessing_parameters.calculate_reynolds_stress    = true;
  postprocessing_parameters.initial_time                 = 0.5;

  TrilinosWrappers::MPI::Vector solution(locally_owned_dofs, mpi_communicator);
  solution(0) = 2.0;
  solution(1) = 0.1;
  solution(2) = 0.0;
  solution(3) = 30;
  solution(4) = 2.5;
  solution(5) = 0.56;
  solution(6) = 0.1;
  solution(7) = 20;


  TrilinosWrappers::MPI::Vector average_solution(locally_owned_dofs,
                                                 mpi_communicator);
  TrilinosWrappers::MPI::Vector stress_solution(locally_owned_dofs,
                                                mpi_communicator);

  // Time info
  const double time_end     = simulation_control_parameters.timeEnd;
  const double initial_time = postprocessing_parameters.initial_time;
  double       time         = simulation_control->get_current_time();
  double       epsilon      = 1e-6;

  while (time < (time_end + epsilon)) // Until time reached end time
    {
      if (time > (initial_time - epsilon)) // Time reached the initial time
        {
          postprocessing_velocities.calculate_average_velocities(
            solution,
            postprocessing_parameters,
            simulation_control->get_current_time(),
            simulation_control->get_time_step(),
            simulation_control->is_output_iteration(),
            locally_owned_dofs,
            mpi_communicator);

          average_solution = postprocessing_velocities.get_average_velocities();
          stress_solution  = postprocessing_velocities.get_reynolds_stress();

          deallog << " Time  : " << time << std::endl;
          deallog << "<u'u'> : " << stress_solution[0] << " "
                  << stress_solution[4] << std::endl;
          deallog << "<v'v'> : " << stress_solution[1] << " "
                  << stress_solution[5] << std::endl;
          deallog << "<w'w'> : " << stress_solution[2] << " "
                  << stress_solution[6] << std::endl;
          deallog << "<u'v'> : " << stress_solution[3] << " "
                  << stress_solution[7] << std::endl;
          deallog << "" << std::endl;
        }

      // New solution values for next step
      solution *= 0.9;

      // Integrate to get the next time
      simulation_control->integrate();

      // Break if the next time from integrate() is the same because
      // time will never get over the time end, but the average velocities
      // at this time is wanted.
      if (abs(time - simulation_control->get_current_time()) < epsilon)
        break;

      time = simulation_control->get_current_time();
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      test(argc, argv);
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

  return 0;
}
