/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

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
  const std::string                                filename,
  const unsigned int                               output_frequency,
  const std::vector<unsigned int>                  boundary_index,
  const double                                     time_step,
  DEM::dem_data_structures<3>::vector_on_boundary &forces_boundary_information,
  DEM::dem_data_structures<3>::vector_on_boundary
    &torques_boundary_information);
#endif
