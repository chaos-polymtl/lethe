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

#include <dem/contact_type.h>
#include <dem/data_containers.h>
#include <dem/dem_container_manager.h>
#include <dem/particle_point_line_contact_info_struct.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

#ifndef locate_local_particles_h
#  define locate_local_particles_h

/**
 * Updates the iterators to particles in local-local contact containers. This is
 * essential since sort_particles_into_subdomains_and_cells() and
 * exchange_ghost_particles() functions change the iterator to particles
 * everytime they are called.
 *
 * @tparam dim Dimensionality of the geometry which contains the particles
 * @param particle_handler
 * @param container_manager The dem container manager that contains
 * particle_container, ghost_adjacent_particles, local_adjacent_particles,
 * particle_wall_pairs_in_contact, particle_floating_wall_pairs_in_contact,
 * particle_floating_mesh_pairs_in_contact, particle_points_in_contact,
 * particle_lines_in_contact
 *
 */

template <int dim>
void
locate_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler,
  DEMContainerManager<dim> &             container_manager);

/**
 * Updates the iterators to local particles in a map of particles
 * (particle_container)
 *
 * @tparam dim Dimensionality of the geometry which contains the particles
 * @param particle_container A map of particles which is used to update
 * the iterators to particles in particle-particle and particle-wall fine search
 * outputs after calling sort particles into cells function
 * @param particle_handler Particle handler to access all the particles in the
 * system
 */

template <int dim>
void
update_particle_container(
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &                                    particle_container,
  const Particles::ParticleHandler<dim> *particle_handler);

/**
 * Updates the iterators to particles in pairs_in_contact
 * (output of particle-object fine search)
 *
 * @tparam dim Dimensionality of the geometry which contains the particles
 * @tparam pairs_structure DEM data structure which contains particle-object pairs relevant information
 * @tparam contact_type Contact type of the contact pairs
 * @param pairs_in_contact Output of particle-object fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim, typename pairs_structure, ContactType contact_type>
void
update_contact_container_iterators(
  pairs_structure &pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container);


#endif /* locate_local_particles_h */
