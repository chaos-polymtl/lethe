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


#ifndef particle_wall_fine_search_h
#define particle_wall_fine_search_h


#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

#include <core/parameters_lagrangian.h>
#include <dem/pw_contact_info_struct.h>

#include <unordered_map>

using namespace dealii;

/**
 * This class is used for fine particle-wall contact search. Fine search
 * is used to find all the particles which are physically in contact with
 * system boundaries and obtain all the required information for calculation
 * of the corresponding particle-wall contact force
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class PWFineSearch
{
public:
  PWFineSearch<dim>();

  /**
   * Iterates over the contact candidates from particle-wall broad search
   * (pw_contact_pair_candidates) to add new contact pairs to the
   * pw_pairs_in_contact container
   *
   * @param pw_contact_pair_candidates The output of particle-wall broad search
   * which shows contact pair candidates
   * @param pw_pairs_in_contact An unordered_map of maps which stores
   * all the particle-wall pairs which are physically in contact, and the
   * contact information in a struct. Note that the size of this unordered map
   * is equal to the number of particles
   */

  void
  particle_wall_fine_search(
    std::unordered_map<
      int,
      std::unordered_map<int,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>,
                                    unsigned int>>> &pw_contact_pair_candidates,
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      &pw_pairs_in_contact);

  /**
   * Iterates over the contact candidates from particle-floating wall broad
   * search (pfw_contact_candidates) to add new contact pairs to the
   * pfw_pairs_in_contact container
   *
   * @param pfw_contact_pair_candidates The output of particle-floating wall
   * broad search which shows contact pair candidates
   * @param floating_wall_properties Properties of floating walls defined in the
   * parameter handler
   * @param simulation_time Simulation time
   * @param pfw_pairs_in_contact An unordered_map of maps which stores
   * all the particle-floating wall pairs which are physically in contact, and
   * the contact information in a struct. Note that the size of this unordered
   * map is equal to the number of particles
   */
  void
  particle_floating_wall_fine_search(
    std::unordered_map<
      int,
      std::unordered_map<int, Particles::ParticleIterator<dim>>>
      &                                               pfw_contact_candidates,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double &                                    simulation_time,
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      &pfw_pairs_in_contact);
};

#endif /* particle_wall_fine_search_h */
