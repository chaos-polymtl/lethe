// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_fine_search_h
#define lethe_particle_particle_fine_search_h

#include <dem/contact_type.h>
#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

using namespace dealii;

/**
 * @brief Selects adjacent_particle_pairs or periodic_adjacent_particle_pairs
 * based on the contact_type
 */
template <int dim, ContactType contact_type>
using adjacent_pairs_for_contact_t = std::conditional_t<
  contact_type == local_periodic_particle_particle ||
    contact_type == ghost_periodic_particle_particle ||
    contact_type == ghost_local_periodic_particle_particle,
  typename DEM::dem_data_structures<dim>::periodic_adjacent_particle_pairs,
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs>;


/**
 * @brief Iterates over a vector of maps (pairs_in_contact) to see if the
 * particles which were in contact in the last time step, are still in contact
 * or not. If they are still in contact it will update the collision info,
 * including tangential displacement, based on new properties of the particle
 * pair, if they are not in contact anymore it will delete the pair from the
 * pairs_in_contact and also its information from pairs_in_contact_info.
 * Then it iterates over the contact candidates from broad search to see if
 * they already exist in the pairs_in_contact or not, if they are not in the
 * pairs_in_contact and have an overlap, the pair will be added to the
 * pairs_in_contact and its contact information will be stored in the
 * corresponding element of the pairs_in_contact_info
 *
 * @param particle_container A container that is used to obtain iterators to
 * particles using their ids
 * @param adjacent_particles A map of maps which stores all the required
 * information for calculation of the contact force of particle pairs. Can
 * either be for non-periodic or periodic contact types
 * @param contact_pair_candidates The output of broad search which shows
 * contact pair candidates
 * @param neighborhood_threshold A value which defines the neighbor particles
 * @param combined_periodic_offsets A vector of tensors of the periodic offsets
 * to change the location of the particles crossing periodic boundaries.
 */
template <int dim, ContactType contact_type>
void
particle_particle_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                                  &particle_container,
  adjacent_pairs_for_contact_t<dim, contact_type> &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                                    &contact_pair_candidates,
  const double                       neighborhood_threshold,
  const std::vector<Tensor<1, dim>> &combined_periodic_offsets = {});

#endif
