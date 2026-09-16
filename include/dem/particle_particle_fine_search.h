// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_fine_search_h
#define lethe_particle_particle_fine_search_h

#include <dem/contact_type.h>
#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

#include <array>

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
 * they already exist in the pairs_in_contact or not. Broad search
 * regenerates this candidate list from cell-neighbor relationships on every
 * call, with no regard to existing contact status, so it is a superset that
 * also includes pairs already present in pairs_in_contact; for those, a
 * cheap map lookup (find()/contains(), never operator[], so a candidate
 * that is not actually in contact never auto-vivifies an empty entry) skips
 * the pair immediately, before recomputing its distance or, for periodic
 * contacts, its nearest periodic image, since the loop above already
 * handled it. This lookup also turns the rare duplicate candidates that
 * broad search can generate for the same pair within a single cycle into a
 * cheap no-op. If a candidate is not already in the pairs_in_contact and
 * has an overlap, the pair will be added to the pairs_in_contact and its
 * contact information will be stored in the corresponding element of the
 * pairs_in_contact_info.
 *
 * @param particle_container A container that is used to obtain iterators to
 * particles using their ids.
 * @param adjacent_particles A map of maps which stores all the required
 * information for calculation of the contact force of particle pairs. Can
 * either be for non-periodic or periodic contact types.
 * @param contact_pair_candidates The output of broad search which shows
 * contact pair candidates.
 * @param neighborhood_threshold A value which defines the neighbor particles.
 * @param periodic_offset_per_direction A tensor whose component d holds the
 * signed period of the domain along direction d (0 if d is not periodic),
 * used to find the nearest periodic image of particles crossing periodic
 * boundaries via the minimum image convention.
 * @param inverse_periodic_offset_per_direction An array whose component d
 * holds the reciprocal of periodic_offset_per_direction's component d (0 if
 * d is not periodic), used to avoid a division when applying the minimum
 * image convention. This is a plain array of independent scale factors, not
 * a Tensor, since it is only ever used component-wise.
 */
template <int dim, ContactType contact_type>
void
particle_particle_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                                  &particle_container,
  adjacent_pairs_for_contact_t<dim, contact_type> &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                                &contact_pair_candidates,
  const double                   neighborhood_threshold,
  const Tensor<1, dim>          &periodic_offset_per_direction         = {},
  const std::array<double, dim> &inverse_periodic_offset_per_direction = {});

#endif
