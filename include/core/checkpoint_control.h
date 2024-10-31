// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_checkpoint_control_h
#define lethe_checkpoint_control_h

#include <core/parameters.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

class CheckpointControl
{
public:
  /**
   * @brief The CheckpointControl is responsible for keeping track of the
   * next checkpoint id to use.
   *
   * @param param Restart struct - Controls writing and reading
   * simulation checkpoints.
   */
  CheckpointControl(Parameters::Restart &param)
    : max_checkpoint_id(2)
    , next_checkpoint_id(max_checkpoint_id)
    , checkpointing_frequency(param.frequency)
    , restart(param.checkpoint)
    , filename(param.filename)
  {}

  /**
   * @brief Return true if a checkpoint should be made at the current time step.
   * @param time_step_number Current time step.
   */
  bool
  is_checkpoint_time_step(const unsigned int time_step_number)
  {
    if (restart && !(time_step_number % checkpointing_frequency))
      {
        increment_checkpoint_id();
        return true;
      }
    return false;
  }

  /**
   * @brief Return the next checkpoint id.
   */
  unsigned int
  get_next_checkpoint_id() const
  {
    return next_checkpoint_id;
  }

  /**
   * @brief Return the prefix for the restart files.
   */
  std::string
  get_filename() const
  {
    return filename;
  }

  /**
   * @brief Serialize the checkpoint controller object to an output archive.
   * @param ar Output archive where the attributes are stored.
   */
  void
  serialize(boost::archive::text_oarchive &ar) const
  {
    ar &next_checkpoint_id;
  }

  /**
   * @brief Deserialize an input archive to the checkpoint controller object.
   * @param ar Input archive where the attributes are stored.
   */
  void
  deserialize(boost::archive::text_iarchive &ar)
  {
    ar &next_checkpoint_id;
  }

private:
  /**
   * @brief Calculate the next checkpoint id.
   */
  void
  increment_checkpoint_id()
  {
    next_checkpoint_id = (next_checkpoint_id + 1) % max_checkpoint_id;
  }

  const unsigned int max_checkpoint_id;
  unsigned int       next_checkpoint_id;
  const unsigned int checkpointing_frequency;
  const bool         restart;
  const std::string  filename;
};
#endif
