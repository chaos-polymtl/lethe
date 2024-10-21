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
   * @brief The CheckpointControl class is responsible for keeping track of the
   * next checkpoint id needed to be use.
   */
  CheckpointControl()
    : next_checkpoint_id(0)
    , max_checkpoint_id(2)
  {}

  /**
   * @brief Increment the next checkpoint id variable by one and apply the modulo.
   */
  void
  increment_checkpoint_id()
  {
    next_checkpoint_id = (next_checkpoint_id + 1) % max_checkpoint_id;
  }

  /**
   * @brief Returns the next checkpoint id.
   */
  unsigned int
  get_next_checkpoint_id()
  {
    return next_checkpoint_id;
  }

  /**
   * @brief Serialize the checkpoint controller object to an output archive.
   * @param ar Output archive where the attributes are stored.
   */
  void
  serialize(boost::archive::text_oarchive &ar, const unsigned int)
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
  unsigned int       next_checkpoint_id;
  const unsigned int max_checkpoint_id;
};
#endif
