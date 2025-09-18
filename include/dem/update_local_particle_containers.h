// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_update_local_particle_containers_h
#define lethe_update_local_particle_containers_h

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/data_containers.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

/**
 * @brief Update the iterators to local particles in a map of particles
 * (particle_container)
 *
 * @tparam dim Dimensionality of the geometry which contains the particles
 *
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
                                        &particle_container,
  const Particles::ParticleHandler<dim> *particle_handler);

/**
 * @brief Update the iterators to particles in pairs_in_contact
 * (output of particle-object fine search)
 *
 * @tparam dim Dimensionality of the geometry which contains the particles
 * @tparam pairs_structure Adjacent particle-object pairs container type.
 * @tparam contact_type Label of contact type to apply proper manipulation of
 * contact removal in containers.
 *
 * @param pairs_in_contact Output of particle-object fine search.
 * @param particle_container Output of update_particle_container function.
 */
template <int dim, typename pairs_structure, ContactType contact_type>
void
update_contact_container_iterators(
  pairs_structure &pairs_in_contact,
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container);


#endif
