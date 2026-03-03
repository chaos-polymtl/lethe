// SPDX-FileCopyrightText: Copyright (c) 2021-2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_output_force_torque_calculation_h
#define lethe_output_force_torque_calculation_h

#include <dem/data_containers.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/table_handler.h>

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
#endif
