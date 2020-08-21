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
#include <dem/pp_contact_info_struct.h>

using namespace dealii;

#ifndef UPDATEGHOSTPPCONTACTCONTAINER_H_
#  define UPDATEGHOSTPPCONTACTCONTAINER_H_

/**
 * Updates the iterators to particles in local_ghost adjacent_particles
 * (output of pp fine search)
 *
 * @param ghost_adjacent_particles Output of particle-particle fine search
 * @param particle_container Output of update_particle_container function
 */

template <int dim>
void
update_ghost_pp_contact_container_iterators(
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container);

#endif /* UPDATEGHOSTPPCONTACTCONTAINER_H_ */
