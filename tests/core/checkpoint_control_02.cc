// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the time and wall_time checkpoint controls: checkpoints at the
 * multiples of the time frequency despite round-off errors in the time, and
 * wall time checkpoints at log iterations only, decided by rank 0.
 */

// Lethe
#include <core/checkpoint_control.h>
#include <core/parameters.h>
#include <core/simulation_control.h>

// Tests (with common definitions)
#include <../tests/tests.h>



void
test()
{
  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

  Parameters::SimulationControl simulation_control_parameters;
  simulation_control_parameters.method =
    Parameters::SimulationControl::TimeSteppingMethod::bdf1;
  simulation_control_parameters.dt       = 0.1;
  simulation_control_parameters.time_end = 100.;

  // Time control: a checkpoint at the first time step that reaches each
  // multiple of 0.25, that is at t = 0.3, 0.5, 0.8 and 1.0. The time of the
  // last one is slightly below 1.0 because of round-off errors.
  {
    SimulationControlTransient simulation_control(
      simulation_control_parameters);
    Parameters::Restart restart_parameters;
    restart_parameters.checkpoint = true;
    restart_parameters.checkpoint_control =
      Parameters::Restart::CheckpointControlType::time;
    restart_parameters.time_frequency = 0.25;
    CheckpointControl checkpoint_control(restart_parameters);

    for (unsigned int time_step = 1; time_step <= 10; ++time_step)
      {
        simulation_control.integrate();
        if (checkpoint_control.is_checkpoint_time_step(simulation_control,
                                                       false,
                                                       mpi_communicator))
          deallog << "Time control: checkpoint at t = "
                  << simulation_control.get_current_time() << std::endl;
      }
  }

  // Wall time control, checked at log iterations (the even time steps): with a
  // wall time frequency that is always elapsed, a checkpoint at every log
  // iteration; with one that is never elapsed, no checkpoint. Only the
  // frequency of rank 0 must matter, so the other ranks have the opposite one.
  simulation_control_parameters.log_frequency = 2;
  for (const double wall_time_frequency : {1e-9, 1e9})
    {
      SimulationControlTransient simulation_control(
        simulation_control_parameters);
      Parameters::Restart restart_parameters;
      restart_parameters.checkpoint = true;
      restart_parameters.checkpoint_control =
        Parameters::Restart::CheckpointControlType::wall_time;
      restart_parameters.wall_time_frequency =
        Utilities::MPI::this_mpi_process(mpi_communicator) == 0 ?
          wall_time_frequency :
          1. / wall_time_frequency;
      CheckpointControl checkpoint_control(restart_parameters);

      for (unsigned int time_step = 1; time_step <= 4; ++time_step)
        {
          simulation_control.integrate();
          if (checkpoint_control.is_checkpoint_time_step(simulation_control,
                                                         false,
                                                         mpi_communicator))
            deallog << "Wall time control (" << wall_time_frequency
                    << " s): checkpoint at time step " << time_step
                    << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      MPILogInitAll                    all;

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

  return 0;
}
