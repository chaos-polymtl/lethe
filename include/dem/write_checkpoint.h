/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <core/pvd_handler.h>
#include <core/serial_solid.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <fstream>
#include <iostream>

using namespace std;

#ifndef write_checkpoint_h
#  define write_checkpoint_h

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
 */

template <int dim>
void
write_checkpoint(
  TimerOutput                                             &computing_timer,
  const DEMSolverParameters<dim>                          &parameters,
  std::shared_ptr<SimulationControl>                      &simulation_control,
  PVDHandler                                              &particles_pvdhandler,
  PVDHandler                                              &grid_pvdhandler,
  parallel::distributed::Triangulation<dim>               &triangulation,
  Particles::ParticleHandler<dim>                         &particle_handler,
  std::shared_ptr<Insertion<dim>>                         &insertion_object,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solid_objects,
  const ConditionalOStream                                &pcout,
  MPI_Comm                                                &mpi_communicator);

#endif /* write_checkpoint_h */
