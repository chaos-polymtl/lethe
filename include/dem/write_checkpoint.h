// SPDX-FileCopyrightText: Copyright (c) 2021, 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_write_checkpoint_h
#define lethe_write_checkpoint_h

#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

/**
 * @brief Write_checkpoint Write a DEM simulation checkpointing to allow for DEM
 * simulation restart.
 *
 * @param computing_timer Dem timer
 * @param parameters Input DEM parameters in the parameter handler file
 * @param simulation_control Simulation control
 * @param particles_pvdhandler PVD handler
 * @param grid_pvdhandler PVD handler for post-processing
 * @param triangulation Triangulation
 * @param particle_handler Particle handler
 * @param insertion_object Insertion object
 * @param solid_objects Vector of solids objects used in DEM simulations
 * @param pcout Printing in parallel
 * @param mpi_communicator
 * @param checkpoint_controller Checkpoint controller
 */
template <int dim, typename PropertiesIndex>
void
write_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim, PropertiesIndex>>        &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_objects,
  const ConditionalOStream                                &pcout,
  MPI_Comm                                                &mpi_communicator,
  const CheckpointControl &checkpoint_controller);

#endif
