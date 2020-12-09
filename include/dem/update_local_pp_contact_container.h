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

#ifndef update_local_particle_particle_contact_container_h
#define update_local_particle_particle_contact_container_h

#include <deal.II/base/config.h>

#include <deal.II/particles/particle_iterator.h>

#include <dem/pp_contact_info_struct.h>

#include <unordered_map>

using namespace dealii;

/**
 * Updates the iterators to particles in local-local adjacent_particles
 * (output of pp fine search)
 *
 * @param local_adjacent_particles Output of particle-particle fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim>
void
update_local_pp_contact_container_iterators(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container);

#endif /* update_local_particle_particle_contact_container_h */
