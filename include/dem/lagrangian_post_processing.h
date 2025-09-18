// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_lagrangian_post_processing_h
#define lethe_lagrangian_post_processing_h


#include <core/pvd_handler.h>

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/distributed/tria.h>

using namespace dealii;

/**
 * This class deals with Lagrangian post-processing. At the moment, calculation
 * of time-averaged particles velocity and granular temperature are added to the
 * class.
 *
 * @todo Add tracer particles
 */

/**
 * @brief Carries out writing the cell data of the domain.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param triangulation Triangulation of the domain.
 * @param grid_pvdhandler PVD handler for grid.
 * @param background_dh Background DoF handler.
 * @param particle_handler Particle handler.
 * @param dem_parameters DEM parameters.
 * @param current_time Simulation time.
 * @param step_number DEM step number.
 * @param mpi_communicator MPI communicator.
 * @param sparse_contacts_object Àdaptive sparse contacts object.
 */
template <int dim, typename PropertiesIndex>
void
write_post_processing_results(
  const parallel::distributed::Triangulation<dim> &triangulation,
  PVDHandler                                      &grid_pvdhandler,
  const DoFHandler<dim>                           &background_dh,
  const Particles::ParticleHandler<dim>           &particle_handler,
  const DEMSolverParameters<dim>                  &dem_parameters,
  const double                                     current_time,
  const unsigned int                               step_number,
  const MPI_Comm                                  &mpi_communicator,
  AdaptiveSparseContacts<dim, PropertiesIndex>    &sparse_contacts_object);

/**
 * @brief Carries out the calculation of the average particles velocity in each local
 * cell. These values are summed up during the post-processing steps (imposed
 * by post-processing frequency), and finally divided by the sampling counter
 * to calculate time-averaged particles velocity distribution.
 *
 * @param triangulation Triangulation.
 * @param particle_handler Particle handler.
 * @param velocity_average_x Average velocity in x-direction.
 * @param velocity_average_y Average velocity in y-direction.
 * @param velocity_average_z Average velocity in z-direction.
 * @param velocity_average_magnitude Average velocity magnitude.
 */
template <int dim, typename PropertiesIndex>
void
calculate_average_particles_velocity(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim>           &particle_handler,
  Vector<double>                                  &velocity_average_x,
  Vector<double>                                  &velocity_average_y,
  Vector<double>                                  &velocity_average_z,
  Vector<double>                                  &velocity_average_magnitude);

/**
 * @brief Carries out the calculation of the granular temperature in each local cell.
 * These values are summed up during the post-processing steps (imposed by
 * post-processing frequency), and finally divided by the sampling counter to
 * calculate time-averaged granular temperature distribution.
 *
 * @param triangulation Triangulation.
 * @param particle_handler Particle handler.
 * @param granular_temperature_average Average granular temperature.
 */
template <int dim, typename PropertiesIndex>
void
calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim>           &particle_handler,
  Vector<double> &granular_temperature_average);

/**
 * @brief Carries out the calculation of average particles velocity for a single
 * cell. It is defined as a separate private function since the
 * calculate_average_granular_temperature function also needs to obtain the
 * average velocity in the cell.
 *
 * @param cell A single cell in the triangulation
 * @param particle_handler
 *
 * @return A tensor which stores the average particles velocity in the input cell
 */
template <int dim, typename PropertiesIndex>
Tensor<1, dim>
calculate_cell_average_particles_velocity(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const Particles::ParticleHandler<dim> &particle_handler);

#endif
