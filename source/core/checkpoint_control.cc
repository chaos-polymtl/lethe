// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/checkpoint_control.h>
#include <core/utilities.h>

#include <cmath>
#include <fstream>

using namespace dealii;

bool
CheckpointControl::is_checkpoint_time_step(
  SimulationControl &simulation_control,
  const bool         is_last_time_step,
  const MPI_Comm    &mpi_communicator)
{
  if (!restart)
    return false;

  bool is_checkpoint = is_last_time_step;
  switch (checkpoint_control)
    {
      case Parameters::Restart::CheckpointControlType::iteration:
        is_checkpoint =
          is_checkpoint || (simulation_control.get_iteration_number() %
                              checkpointing_frequency ==
                            0);
        break;
      case Parameters::Restart::CheckpointControlType::time:
        {
          // A checkpoint is written when the time crosses a multiple of the
          // time frequency. The tolerance, the same as
          // SimulationControlTransient::output_time_tolerance, prevents
          // round-off errors in the time from delaying or skipping a
          // checkpoint.
          constexpr double time_tolerance = 1e-6;
          const double     margin         = time_tolerance * time_frequency;
          is_checkpoint =
            is_checkpoint ||
            std::floor((simulation_control.get_current_time() + margin) /
                       time_frequency) >
              std::floor((simulation_control.get_previous_time() + margin) /
                         time_frequency);
          break;
        }
      case Parameters::Restart::CheckpointControlType::wall_time:
        // The wall time differs between ranks, so rank 0 decides, at log
        // iterations only.
        if (simulation_control.is_verbose_iteration())
          is_checkpoint =
            Utilities::MPI::broadcast(mpi_communicator,
                                      wall_time_timer.wall_time() >=
                                        wall_time_frequency,
                                      0) ||
            is_checkpoint;
        break;
    }

  if (!is_checkpoint)
    return false;

  increment_checkpoint_id();
  wall_time_timer.restart();
  return true;
}

void
CheckpointControl::write(const std::string &prefix,
                         const MPI_Comm    &mpi_communicator) const
{
  // This file has the same name for every checkpoint id. It gives the id of
  // the last checkpoint, from which a restart reads.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::ofstream                 output(prefix + ".checkpoint_controller");
      boost::archive::text_oarchive archive(output, boost::archive::no_header);
      serialize(archive);
    }
}

void
CheckpointControl::read(const std::string &prefix,
                        const MPI_Comm    &mpi_communicator)
{
  next_checkpoint_id =
    read_checkpoint_id(prefix,
                       "checkpoint file of the checkpoint controller, given by "
                       "'subsection restart' - 'set filename'",
                       mpi_communicator);
}

unsigned int
CheckpointControl::read_checkpoint_id(const std::string &prefix,
                                      const std::string &file_description,
                                      const MPI_Comm    &mpi_communicator)
{
  const std::string filename = prefix + ".checkpoint_controller";

  // Checked on every rank, so that a missing file throws on every rank instead
  // of leaving the other ranks waiting in the broadcast.
  check_file_exists(filename, file_description);

  unsigned int checkpoint_id = 0;
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::ifstream                 input(filename);
      boost::archive::text_iarchive archive(input, boost::archive::no_header);
      archive >> checkpoint_id;
    }
  return Utilities::MPI::broadcast(mpi_communicator, checkpoint_id, 0);
}
