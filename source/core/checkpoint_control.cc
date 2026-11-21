// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/checkpoint_control.h>
#include <core/utilities.h>

#include <fstream>

using namespace dealii;

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
