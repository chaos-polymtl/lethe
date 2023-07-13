/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <core/dem_properties.h>
#include <core/serial_solid.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

#ifndef find_contact_detection_step_h
#  define find_contact_detection_step_h

/**
 * @brief Carries out finding steps for dynamic contact search for particle-particle
 * contacts
 *
 * @param particle_handler
 * @param dt DEM time step
 * @param smallest_contact_search_criterion A criterion for finding
 * dynamic contact search steps. This value is defined as the minimum of
 * particle-particle and particle-wall displacement threshold values
 * @param mpi_communicator
 * @param sorting_in_subdomains_step True if it is insertion, load-
 * balance or contact detection step
 * @param displacement Displacement of particles since last sorting step
 * @return Returns true if the maximum cumulative
 * displacement of particles exceeds the threshold and false otherwise
 *
 */

template <int dim>
bool
find_particle_contact_detection_step(
  Particles::ParticleHandler<dim> &particle_handler,
  const double                     dt,
  const double                     smallest_contact_search_criterion,
  MPI_Comm &                       mpi_communicator,
  bool                             sorting_in_subdomains_step,
  std::vector<double> &            displacement);

/**
 * @brief Carries out finding steps for dynamic contact search in particle-floating
 * mesh contacts
 *
 * @param smallest_contact_search_criterion A criterion which defines the maximal displacement that a solid face may have displaced
 * @param solids All solid objects used in the simulation
 * @return Returns true if the maximum cumulative
 * displacement of particles exceeds the threshold and false otherwise
 *
 */

template <int dim>
bool
find_floating_mesh_mapping_step(
  const double smallest_contact_search_criterion,
  std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids);

#endif
