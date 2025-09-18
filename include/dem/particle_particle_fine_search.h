// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_fine_search_h
#define lethe_particle_particle_fine_search_h

#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

using namespace dealii;

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
 * information for calculation of the contact force of particle pairs
 * @param contact_pair_candidates The output of broad search which shows
 * contact pair candidates
 * @param neighborhood_threshold A value which defines the neighbor particles
 * @param periodic_offset A tensor of the periodic offset to change the
 * particle location of the particles on the periodic boundary 1 side,
 * the tensor as 0.0 values by default
 */
template <int dim>
void
particle_particle_fine_search(
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                      &contact_pair_candidates,
  const double         neighborhood_threshold,
  const Tensor<1, dim> periodic_offset = Tensor<1, dim>());

#endif
