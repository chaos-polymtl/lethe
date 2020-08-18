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

using namespace dealii;

#ifndef UPDATELOCALPPCONTACTCONTAINER_H_
#define UPDATELOCALPPCONTACTCONTAINER_H_

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
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
    &local_adjacent_particles,
  const std::unordered_map<int, Particles::ParticleIterator<dim>>
    &particle_container)
{
  for (auto adjacent_particles_iterator = local_adjacent_particles.begin();
       adjacent_particles_iterator != local_adjacent_particles.end();
       ++adjacent_particles_iterator)
    {
      int  particle_one_id          = adjacent_particles_iterator->first;
      auto pairs_in_contant_content = &adjacent_particles_iterator->second;
      for (auto pp_map_iterator = pairs_in_contant_content->begin();
           pp_map_iterator != pairs_in_contant_content->end();
           ++pp_map_iterator)
        {
          int particle_two_id = pp_map_iterator->first;

          pp_map_iterator->second.particle_one =
            particle_container.at(particle_one_id);
          pp_map_iterator->second.particle_two =
            particle_container.at(particle_two_id);
        }
    }
}

#endif /* UPDATELOCALPPCONTACTCONTAINER_H_ */
