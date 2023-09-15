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
#include <dem/particle_particle_contact_info.h>
#include <dem/particle_wall_contact_info.h>

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

template <int dim,
          typename pairs_structure,
          typename candidates_structure,
          ContactType contact_type>
void
update_fine_search_candidates(pairs_structure      &pairs_in_contact,
                              candidates_structure &contact_candidates);

#endif /* localize_contacts_h */
