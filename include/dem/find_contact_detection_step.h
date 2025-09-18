// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_find_contact_detection_step_h
#define lethe_find_contact_detection_step_h

#include <core/serial_solid.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

/**
 * @brief Find steps for dynamic contact search for particle-particle contacts.

 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle_handler
 * @param dt DEM time step
 * @param smallest_contact_search_criterion A criterion for finding
 * dynamic contact search steps. This value is defined as the minimum of
 * particle-particle and particle-wall displacement threshold values
 * @param mpi_communicator
 * @param displacement Displacement of particles since last sorting step
 * @param parallel_update Update the identification of the contact detection
 * step in parallel. If this parameter is set to false, the distance will be
 * calculated but the
 * logical OR statement won't be called and a false value will be returned. In
 * essence, this will only update the displacement.
 *
 * @return Returns true if the maximum cumulative displacement of particles
 * exceeds the threshold and false otherwise

 */
template <int dim, typename PropertiesIndex>
void
find_particle_contact_detection_step(
  Particles::ParticleHandler<dim> &particle_handler,
  const double                     dt,
  const double                     smallest_contact_search_criterion,
  MPI_Comm                        &mpi_communicator,
  std::vector<double>             &displacement,
  const bool                       parallel_update = true);

/**
 * @brief Find steps for dynamic contact search in particle-floating
 * mesh contacts
 *
 * @param smallest_contact_search_criterion A criterion which defines the maximal displacement that a solid face may have displaced
 * @param solids All solid objects used in the simulation
 * @return Returns true if the maximum cumulative
 * displacement of particles exceeds the threshold and false otherwise
 */
template <int dim>
void
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids);

#endif
