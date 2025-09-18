// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_read_checkpoint_h
#define lethe_read_checkpoint_h

#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

/**
 * @brief Read a DEM simulation checkpoint, allowing the simulation to restart
 * from where it stopped.
 *
 * @param computing_timer Dem timer
 * @param parameters Input DEM parameters in the parameter handler file
 * @param simulation_control Simulation control
 * @param particles_pvdhandler PVD handler
 * @param grid_pvdhandler PVD handler for post-processing
 * @param triangulation Triangulation
 * @param particle_handler Particle handler
 * @param insertion_object Shared pointer of Insertion type.
 * @param solid_surfaces Vector of solids surfaces used in DEM simulations
 * @param checkpoint_controller Checkpoint controller
 */
template <int dim, typename PropertiesIndex>
void
read_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim, PropertiesIndex>>        &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_surfaces,
  CheckpointControl &checkpoint_controller);

/**
 * @brief Asserts that a file exists. If the file does not exist, an exception is
 * @param[in] f Object associated with the file to check. The file
 * should be opened before passing it to this function.
 * @param[in] filename Name of the file that is not found.
 *
 */
inline void
assert_restart_file_exists(std::ifstream &f, const std::string &filename)
{
  AssertThrow(f.is_open(),
              ExcMessage(
                std::string("You are trying to restart a previous computation, "
                            "but the restart file <") +
                filename + "> does not appear to exist!"));
}
#endif
