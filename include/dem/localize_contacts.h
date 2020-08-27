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
#include <dem/pw_contact_info_struct.h>

using namespace std;

#ifndef LOCALIZECONTACTS_H_
#  define LOCALIZECONTACTS_H_

/**
 * Manages clearing the contact containers when particles are exchanged
 * between processors. If the adjacent pair does not exist in the output of the
 * new (at this step) broad search, it is deleted from the adjacent pair, since
 * it means that the contact is being handled on another processor. If the pair
 * exists in the output of the new broad search, it is deleted from the output
 * of the broad search. This process is performed for local-local
 * particle-particle pairs, local-ghost particle-particle pairs and
 * particle-wall pairs.
 *
 * @param local_adjacent_particles Local-local adjacent particle pairs
 * @param ghost_adjacent_particles Local-ghost adjacent particle pairs
 * @param pw_pairs_in_contact Particle-wall contact pairs
 * @param local_contact_pair_candidates Outputs of local-local particle-particle
 * broad search
 * @param ghost_contact_pair_candidates Outputs of local-ghost particle-particle
 * broad search
 * @param pw_contact_candidates Outputs of particle-wall broad search
 *
 */

template <int dim>
void
localize_contacts(
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *local_adjacent_particles,
  std::unordered_map<int, std::unordered_map<int, pp_contact_info_struct<dim>>>
    *ghost_adjacent_particles,
  std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
    *                                        pw_pairs_in_contact,
  std::unordered_map<int, std::vector<int>> &local_contact_pair_candidates,
  std::unordered_map<int, std::vector<int>> &ghost_contact_pair_candidates,
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    &pw_contact_candidates);

#endif /* LOCALIZECONTACTS_H_ */
