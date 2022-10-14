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
 */

#include <dem/contact_type.h>
#include <dem/data_containers.h>

using namespace dealii;

#ifndef update_fine_search_contacts_h
#  define update_fine_search_contacts_h

/**
 *
 *
 * @param
 * @param
 */

template <int dim, typename pairs_structure, typename candidates_structure>
void
update_fine_search_candidates(pairs_structure &     pairs_in_contact,
                              candidates_structure &contact_candidates,
                              const ContactType     contact_type);

template <int dim>
void
update_particle_fine_search_candidates(
  typename dem_data_containers::dem_data_structures<
    dim>::adjacent_particle_pairs &adjacent_particles,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_particle_candidates &contact_pair_candidates,
  const ContactType                     contact_type);


template <int dim, typename pairs_structure, typename candidates_structure>
void
update_wall_fine_search_candidates(
  pairs_structure &     particle_wall_pairs_in_contact,
  candidates_structure &particle_wall_contact_candidates);


template <int dim>
void
update_mesh_fine_search_candidates(
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
  typename dem_data_containers::dem_data_structures<
    dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates);


#endif /* update_fine_search_contacts_h */