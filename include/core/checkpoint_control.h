// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_checkpoint_control_h
#define lethe_checkpoint_control_h

#include <core/parameters.h>
#include <core/simulation_control.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/timer.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <string>

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
  CheckpointControl(const Parameters::Restart &param)
    : max_checkpoint_id(2)
    , next_checkpoint_id(0)
    , checkpointing_frequency(param.frequency)
    , restart(param.checkpoint)
    , filename(param.filename)
    , checkpoint_control(param.checkpoint_control)
    , time_frequency(param.time_frequency)
    , wall_time_frequency(param.wall_time_frequency)
  {}

  /**
   * @brief Return true if a checkpoint should be made at the current time step:
   * at the last time step, and according to the checkpoint control at each
   * checkpointing frequency time step, each time frequency of simulated time,
   * or at the first log iteration after each wall time frequency.
   * @param simulation_control Simulation control of the simulation.
   * @param is_last_time_step Whether the current time step is the last one.
   * @param mpi_communicator MPI communicator of the simulation.
   * @return Whether a checkpoint is written at this time step.
   */
  bool
  is_checkpoint_time_step(SimulationControl &simulation_control,
                          const bool         is_last_time_step,
                          const MPI_Comm    &mpi_communicator);

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
   * @brief Write the id of the last checkpoint in the file
   * <prefix>.checkpoint_controller, from which a restart reads. Only rank 0
   * writes it. It is called once every file of the checkpoint is written.
   *
   * @param[in] prefix Prefix of the checkpoint files, without checkpoint slot.
   * @param[in] mpi_communicator MPI communicator of the simulation.
   */
  void
  write(const std::string &prefix, const MPI_Comm &mpi_communicator) const;

  /**
   * @brief Read the id of the last checkpoint from the file
   * <prefix>.checkpoint_controller. It is called before the checkpoint is read.
   *
   * @param[in] prefix Prefix of the checkpoint files, without checkpoint slot.
   * @param[in] mpi_communicator MPI communicator of the simulation.
   */
  void
  read(const std::string &prefix, const MPI_Comm &mpi_communicator);

  /**
   * @brief Return the id of the last checkpoint recorded in the file
   * <prefix>.checkpoint_controller. Only rank 0 reads the file, and the id is
   * broadcast to the other ranks.
   *
   * @param[in] prefix Prefix of the checkpoint files, without checkpoint slot.
   * @param[in] file_description Description of the file, and of the parameter
   * that gives it, reported if the file cannot be read.
   * @param[in] mpi_communicator MPI communicator of the simulation.
   *
   * @return Id of the last checkpoint.
   */
  static unsigned int
  read_checkpoint_id(const std::string &prefix,
                     const std::string &file_description,
                     const MPI_Comm    &mpi_communicator);

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

  const Parameters::Restart::CheckpointControlType checkpoint_control;
  const double                                     time_frequency;
  const double                                     wall_time_frequency;

  // Wall time since the last checkpoint, for the wall_time checkpoint control
  Timer wall_time_timer;
};
#endif
