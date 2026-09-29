// SPDX-FileCopyrightText: Copyright (c) 2021-2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_output_force_torque_calculation_h
#define lethe_output_force_torque_calculation_h

#include <dem/data_containers.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/tensor.h>

#include <string>

/**
 * @brief write_forces_torques_output_locally
 * Writes the results of force and torque calculations in the terminal at the
 * frequency requested in the prm.
 */
void
write_forces_torques_output_locally(
  std::map<unsigned int, Tensor<1, 3>> force_on_walls,
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls);

/**
 * @brief write_forces_torques_output_results
 * Writes the results of force and torque calculations in a file, and it depends
 * on the verbosity, in the terminal
 */
void
write_forces_torques_output_results(
  const std::string                               &filename,
  const unsigned int                               output_frequency,
  const std::vector<unsigned int>                 &boundary_index,
  const double                                     time_step,
  DEM::dem_data_structures<3>::vector_on_boundary &forces_boundary_information,
  DEM::dem_data_structures<3>::vector_on_boundary
    &torques_boundary_information);

/**
 * @brief Prepare the file in which the force and torque exerted by the
 * particles on a solid object are written. For a new simulation, the file is
 * created, or overwritten, and only contains the header. For a restarted
 * simulation, the rows of the existing file whose time is later than the
 * checkpoint time are removed, so that the rows appended afterward are neither
 * duplicated nor out of order. If the file does not exist, it is created as for
 * a new simulation. This function must be called by a single process.
 *
 * @param[in] filename Name of the file.
 * @param[in] restart Whether the simulation is restarted from a checkpoint.
 * @param[in] restart_time Time of the checkpoint the simulation is restarted
 * from. It is not used for a new simulation.
 * @param[in] time_tolerance Tolerance used to compare the time of the rows with
 * the checkpoint time. It is not used for a new simulation.
 */
void
initialize_solid_forces_torques_file(const std::string &filename,
                                     const bool         restart,
                                     const double       restart_time,
                                     const double       time_tolerance);

/**
 * @brief Append a row with the force and torque exerted by the particles on a
 * solid object to the file prepared by initialize_solid_forces_torques_file().
 * The row contains the time followed by the components of the force and of
 * the torque. This function must be called by a single process.
 *
 * @param[in] filename Name of the file.
 * @param[in] time Time of the force and torque.
 * @param[in] force Force exerted by the particles on the solid object.
 * @param[in] torque Torque exerted by the particles on the solid object about
 * its center of rotation.
 */
void
append_solid_forces_torques_to_file(const std::string  &filename,
                                    const double        time,
                                    const Tensor<1, 3> &force,
                                    const Tensor<1, 3> &torque);
#endif
