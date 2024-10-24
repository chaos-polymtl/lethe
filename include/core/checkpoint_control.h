// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_checkpoint_control_h
#define lethe_checkpoint_control_h

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

class CheckpointControl
{
public:
  /**
   * @brief The CheckpointControl is responsible for keeping track of the
   * next checkpoint id to use.
   *
   * @param checkpoint_freq Checkpoint frequency.
   * @param restart_filename Filename used for the prefix during checkpointing.
   */
  CheckpointControl(const unsigned int checkpoint_freq,
                    const std::string &restart_filename)
    : max_checkpoint_id(2)
    , next_checkpoint_id(max_checkpoint_id)
    , checkpointing_frequency(checkpoint_freq)
    , filename(restart_filename)
  {}

  /**
   * @brief Return true if a checkpoint should be made at the current time step.
   * @param time_step_number Current time step.
   */
  bool
  is_checkpoint_time_step(const unsigned int time_step_number)
  {
    return (time_step_number % checkpointing_frequency) == 0;
  }

  /**
   * @brief Calculate the next checkpoint id.
   */
  void
  increment_checkpoint_id()
  {
    next_checkpoint_id = (next_checkpoint_id + 1) % max_checkpoint_id;
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
   * @brief Return de prefix for the restart files.
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
  serialize(boost::archive::text_oarchive &ar, const unsigned int) const
  {
    ar &next_checkpoint_id;
  }

  /**
   * @brief Deserialize an input archive to the checkpoint controller object.
   * @param ar Input archive where the attributes are stored.
   */
  void
  deserialize(boost::archive::text_iarchive &ar, const unsigned int)
  {
    ar &next_checkpoint_id;
  }

private:
  const unsigned int max_checkpoint_id;
  unsigned int       next_checkpoint_id;
  const unsigned int checkpointing_frequency;
  const std::string  filename;
};
#endif
