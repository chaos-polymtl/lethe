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
#include <dem/particle_particle_contact_info_struct.h>
#include <dem/particle_wall_contact_info_struct.h>

using namespace std;

#ifndef localize_contacts_h
#  define localize_contacts_h

/**
 * Manages to call update_fine_search_candidates() with contact data containers
 * to remove contact repetitions and to add new contact pairs to the contact
 * containers when particles are exchanged between processors. Repeats the
 * update for local particle-particle contacts, ghost particle-particle,
 * particle-wall contacts, particle-floating wall contacts and particle-floating
 * mesh contacts.
 *
 * @param local_adjacent_particles Local-local adjacent particle pairs
 * @param ghost_adjacent_particles Local-ghost adjacent particle pairs
 * @param particle_wall_pairs_in_contact Particle-wall contact pairs
 * @param particle_floating_wall_in_contact Particle-floating wall contact pairs
 * @param particle_floating_mesh_in_contact Particle-floating mesh contact pairs
 * @param local_contact_pair_candidates Local-local particle-particle contact candidates
 * @param ghost_contact_pair_candidates Local-ghost particle-particle contact candidates
 * @param particle_wall_contact_candidates Particle-wall contact candidates
 * @param particle_floating_wall_contact_candidates Particle-floating wall contact candidates
 * @param particle_floating_mesh_contact_candidates Particle-floating mesh contact candidates
 *
 */

template <int dim>
void
localize_contacts(
  typename dem_data_containers::dem_data_structures<
    dim>::adjacent_particle_pairs &local_adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    dim>::adjacent_particle_pairs &ghost_adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_in_contact &particle_wall_pairs_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_in_contact &particle_floating_wall_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_particle_candidates &local_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_particle_candidates &ghost_contact_pair_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_wall_candidates &particle_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_wall_candidates
    &particle_floating_wall_contact_candidates,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates);

/**
 * Manages removing repetitions and adding new contact pairs to the contact
 * containers when particles are exchanged between processors. If the contact
 * pair (in adjacent particles containers) does not exist in the output of the
 * new (current step) broad search, it is removed from the contact pair
 * (adjacent containers), since it means that the contact is being handled by
 * another processor. If the pair exists in the output of the new broad search,
 * it is removed from the output of the broad search, as the contact is already
 * being processed. This process is performed for many types of particle-object
 * pairs : local-local particle-particle pairs, local-ghost particle-particle
 * pairs and particle-wall pairs, particle-floating wall pairs and
 * particle-floating mesh contacts pairs. Note : contact_type is important for
 * local or ghost particle-particle contact only
 *
 * @param pairs_in_contact adjacent particle-object pairs from fine search at last time step
 * @param contact_candidates particle-object contact pairs from board search at current time step
 * @param contact_type label of contact type to apply proper manipulation of contact removal in containers
 */

template <int dim, typename pairs_structure, typename candidates_structure>
void
update_fine_search_candidates(pairs_structure &     pairs_in_contact,
                              candidates_structure &contact_candidates,
                              const ContactType     contact_type);

#endif /* localize_contacts_h */
