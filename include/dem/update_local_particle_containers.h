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
#include <dem/particle_point_line_contact_info_struct.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

#ifndef locate_local_particles_h
#  define locate_local_particles_h

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
 * @param clear_contact_structures If true, the contact structures will be cleared
 */

template <int dim, typename pairs_structure, ContactType contact_type>
void
update_contact_container_iterators(
  pairs_structure &pairs_in_contact,
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &        particle_container,
  const bool clear_contact_structures = false);


#endif /* locate_local_particles_h */
