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
#include <dem/data_containers.h>
#include <dem/particle_wall_contact_info_struct.h>

using namespace dealii;

#ifndef update_particle_wall_contact_container_h
#  define update_particle_wall_contact_container_h

/**
 * Updates the iterators to particles in particle_wall_contact_container (output
 * of particle-wall fine search)
 *
 * @param particle_wall_pairs_in_contact Output of particle-wall fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim>
void
update_particle_wall_contact_container_iterators(
  std::unordered_map<
    types::particle_index,
    std::map<unsigned int, particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container);

/**
 * Updates the iterators to particles in particle_moving_wall_contact_container
 *
 * @param particle_moving_wall_pairs_in_contact Output of particle-moving wall fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim>
void
update_particle_moving_wall_contact_container_iterators(
  std::unordered_map<types::global_cell_index,
                     std::unordered_map<types::particle_index,
                                        particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<types::particle_index, Particles::ParticleIterator<dim>>
    &particle_container);


#endif /* update_particle_wall_contact_container_h */
