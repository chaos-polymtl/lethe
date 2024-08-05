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

#ifndef lethe_read_checkpoint_h
#define lethe_read_checkpoint_h

#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

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
 */
template <int dim>
void
read_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim>>                         &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_surfaces);

#endif
