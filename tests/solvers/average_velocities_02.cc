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
 * @brief This code tests averaging values in time with Trilinos vectors with
 * adaptative time step scaling.
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
  simulation_control_parameters.dt                           = 0.1;
  simulation_control_parameters.timeEnd                      = 1.0;
  simulation_control_parameters.output_frequency             = 1;
  simulation_control_parameters.adapt                        = true;
  simulation_control_parameters.adaptative_time_step_scaling = 0.95;

  // Variables for AverageVelocities
  AverageVelocities<3, TrilinosWrappers::MPI::Vector, IndexSet> average;

  auto simulation_control =
    std::make_shared<SimulationControlTransient>(simulation_control_parameters);

  IndexSet locally_owned_dofs;
  locally_owned_dofs.add_range(0, 3);

  Parameters::PostProcessing postprocessing_parameters;
  postprocessing_parameters.calculate_average_velocities = true;
  postprocessing_parameters.initial_time                 = 0.5;

  TrilinosWrappers::MPI::Vector solution(locally_owned_dofs, mpi_communicator);
  solution(0) = 0.0;
  solution(1) = 2.5;
  solution(2) = 10;

  TrilinosWrappers::MPI::Vector average_solution(locally_owned_dofs,
                                                 mpi_communicator);

  // Time and output info
  const double time_end     = simulation_control_parameters.timeEnd;
  const double initial_time = postprocessing_parameters.initial_time;
  double       time         = simulation_control->get_current_time();
  double       dt           = 0.0;

  const double epsilon = 1e-6;

  while (time < (time_end + epsilon)) // Until time reached end time
    {
      if (time > (initial_time - epsilon)) // Time reached the initial time
        {
          average.calculate_average_velocities(solution,
                                               postprocessing_parameters,
                                               time,
                                               dt,
                                               locally_owned_dofs,
                                               mpi_communicator);

          average_solution = average.get_average_velocities();

          deallog << " Time :             " << time << std::endl;
          deallog << " Time step :        " << dt << std::endl;
          deallog << " Average solution : " << average_solution[0] << " "
                  << average_solution[1] << " " << average_solution[2]
                  << std::endl;
          deallog << "" << std::endl;
        }

      // new solution values for next step
      solution *= 0.9;

      // Integrate to get the next time
      simulation_control->integrate();

      // Break if the next time from integrate() is the same because
      // time will never get over the time end, but the average velocities
      // at this time is wanted.
      if (abs(time - simulation_control->get_current_time()) < epsilon)
        break;

      dt   = simulation_control->get_time_step();
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
