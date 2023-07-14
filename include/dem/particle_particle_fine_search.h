/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/particle_particle_contact_info.h>

#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>

#include <boost/range/adaptor/map.hpp>

template <int dim>
class DEMContactManager;

using namespace dealii;

#ifndef particle_particle_fine_search_h
#  define particle_particle_fine_search_h

/**
 * This class is used for local-local and local-ghost fine particle-particle
 * contact search. Fine search is used to find all the particle pairs which are
 * physically in contact and obtain all the required information for calculation
 * of the contact force
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleParticleFineSearch
{
public:
  ParticleParticleFineSearch<dim>();

  /**
   * Iterates over a vector of maps (pairs_in_contact) to see if the particles
   * which were in contact in the last time step, are still in contact or not.
   * If they are still in contact it will update the collision info, including
   * tangential overlap, based on new properties of the particle pair, if they
   * are not in contact anymore it will delete the pair from the
   * pairs_in_contact and also its information from pairs_in_contact_info.
   * Then it iterates over the contact candidates from broad search to see if
   * they already exist in the pairs_in_contact or not, if they are not in the
   * pairs_in_contact and they have overlap, the pair will be added to the
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

  void
  particle_particle_fine_search(
    typename DEM::dem_data_structures<dim>::particle_index_iterator_map
      &particle_container,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &adjacent_particles,
    const typename DEM::dem_data_structures<dim>::particle_particle_candidates
      &                  contact_pair_candidates,
    const double         neighborhood_threshold,
    const Tensor<1, dim> periodic_offset = Tensor<1, dim>());
};

#endif /* particle_particle_fine_search_h */
