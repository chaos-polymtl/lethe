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

#include <dem/particle_particle_contact_info_struct.h>
#include <dem/particle_wall_contact_info_struct.h>

using namespace std;

#ifndef localize_contacts_h
#  define localize_contacts_h

/**
 * Manages removing repetitions and adding new contact pairs to the contact
 * containers when particles are exchanged between processors. If the contact
 * pair (in adjacent particles containers) does not exist in the output of the
 * new (current step) broad search, it is removed from the contact pair
 * (adjacent containers), since it means that the contact is being handled by
 * another processor. If the pair exists in the output of the new broad search,
 * it is removed from the output of the broad search, as the contact is already
 * being processed. This process is performed for local-local particle-particle
 * pairs, local-ghost particle-particle pairs and particle-wall pairs.
 *
 * @param local_adjacent_particles Local-local adjacent particle pairs
 * @param ghost_adjacent_particles Local-ghost adjacent particle pairs
 * @param particle_wall_pairs_in_contact Particle-wall contact pairs
 * @param pfw_pairs_in_contact Particle-floating wall contact pairs
 * @param particle_moving_mesh_in_contact Particle-moving mesh contact pairs
 * @param local_contact_pair_candidates Local-local particle-particle contact candidates
 * @param ghost_contact_pair_candidates Local-ghost particle-particle contact candidates
 * @param particle_wall_contact_candidates Particle-wall contact candidates
 * @param pfw_contact_candidates Particle-floating wall contact candidates
 * @param particle_moving_mesh_contact_candidates Particle-moving mesh contact candidates
 *
 */

template <int dim>
void
localize_contacts(
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &local_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<types::particle_index,
                       particle_particle_contact_info_struct<dim>>>
    &ghost_adjacent_particles,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &particle_wall_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<types::particle_index, particle_wall_contact_info_struct<dim>>>
    &pfw_pairs_in_contact,
  std::unordered_map<
    types::particle_index,
    std::map<unsigned int, particle_wall_contact_info_struct<dim>>>
    &particle_moving_mesh_in_contact,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &local_contact_pair_candidates,
  std::unordered_map<types::particle_index, std::vector<types::particle_index>>
    &ghost_contact_pair_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<unsigned int,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  unsigned int,
                                  unsigned int>>>
    &particle_wall_contact_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<unsigned int, Particles::ParticleIterator<dim>>>
    &pfw_contact_candidates,
  std::unordered_map<
    types::particle_index,
    std::unordered_map<
      unsigned int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    &particle_moving_mesh_contact_candidates);

#endif /* localize_contacts_h */
