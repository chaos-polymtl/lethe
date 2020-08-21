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

#include <dem/update_ghost_iterator_pp_contact_container.h>
#include <dem/update_ghost_particle_container.h>

using namespace dealii;

#ifndef LOCATEGHOSTPARTICLES_H_
#  define LOCATEGHOSTPARTICLES_H_

/**
 * Updates the iterators to particles in local-ghost contact containers. This is
 * essential since sort_particles_into_subdomains_and_cells() and
 * exchange_ghost_particles() functions change the iterator to particles
 * everytime they are called.
 *
 * @param particle_handler
 * @param particle_container A container that contains the updated iterators to
 * all local and ghost particles
 * @param ghost_adjacent_particles Container that contains all the contact
 * information of adjacent local-ghost particles
 * @param local_adjacent_particles Container that contains all the contact
 * information of adjacent local-local particles
 *
 */

template <int dim>
void
locate_ghost_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler,
  std::unordered_map<int, Particles::ParticleIterator<dim>>
    &ghost_particle_container,
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
    &ghost_adjacent_particles);

#endif /* LOCATEGHOSTPARTICLES_H_ */
