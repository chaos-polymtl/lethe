// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the checkpoint slots of CheckpointControl: a checkpoint every
 * checkpointing frequency time steps and at the last time step, alternating
 * between the slots 1 and 0, and the slot of the first checkpoint after a
 * restart.
 */

// Lethe
#include <core/checkpoint_control.h>
#include <core/parameters.h>

// Boost
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <sstream>


void
test()
{
  Parameters::Restart restart_parameters;
  restart_parameters.checkpoint = true;
  restart_parameters.frequency  = 2;

  // Checkpoints every 2 time steps and at the last one, the fifth, so that the
  // third checkpoint is written in the slot of the first one.
  CheckpointControl  checkpoint_control(restart_parameters);
  const unsigned int n_time_steps = 5;
  for (unsigned int time_step = 1; time_step <= n_time_steps; ++time_step)
    if (checkpoint_control.is_checkpoint_time_step(time_step,
                                                   time_step == n_time_steps))
      deallog << "Time step " << time_step << ": checkpoint in slot "
              << checkpoint_control.get_next_checkpoint_id() << std::endl;

  // A restart reads the slot of the last checkpoint, and its first checkpoint
  // is written in the other slot.
  std::stringstream stream;
  {
    boost::archive::text_oarchive archive(stream, boost::archive::no_header);
    checkpoint_control.serialize(archive);
  }
  CheckpointControl restarted_checkpoint_control(restart_parameters);
  {
    boost::archive::text_iarchive archive(stream, boost::archive::no_header);
    restarted_checkpoint_control.deserialize(archive);
  }
  deallog << "Restart from slot "
          << restarted_checkpoint_control.get_next_checkpoint_id() << std::endl;
  if (restarted_checkpoint_control.is_checkpoint_time_step(6, false))
    deallog << "Time step 6: checkpoint in slot "
            << restarted_checkpoint_control.get_next_checkpoint_id()
            << std::endl;

  // No checkpoint is written if checkpointing is disabled.
  restart_parameters.checkpoint = false;
  CheckpointControl disabled_checkpoint_control(restart_parameters);
  deallog << "Checkpoint with checkpointing disabled: "
          << disabled_checkpoint_control.is_checkpoint_time_step(2, true)
          << std::endl;
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

  return 0;
}
