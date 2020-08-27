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

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_properties.h>
#include <dem/pw_contact_info_struct.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef PWFINESEARCH_H_
#  define PWFINESEARCH_H_

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
   * Iterates over a vector of maps (pw_pairs_in_contact) to see if the
   * particle-wall pairs which were in contact in the last time step, are still
   * in contact or not. If they are still in contact it will update the
   * collision info, including tangential overlap, based on new properties of
   * the particle and wall, if they are not in contact anymore it will delete
   * the pair from the pw_pairs_in_contact. Then it iterates over the contact
   * candidates from particle-wall broad search to see if they already exist in
   * the pw_pairs_in_contact or not, if they are not in the pw_pairs_in_contact
   * and they have overlap, the pair will be added to the pw_pairs_in_contact
   * and its contact information will be stored
   *
   * @param pw_contact_pair_candidates The output of particle-wall broad search
   * which shows contact pair candidates
   * @param pw_pairs_in_contact A vector of maps which stores all the
   * particle-wall pairs which are physically in contact, and the contact
   * information in a struct. Note that the size of this vector is equal to the
   * number of particles while the key of map (each element of the vector) is
   * the boundary id
   */

  void
  pw_Fine_Search(
    std::unordered_map<
      int,
      std::unordered_map<int,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>>>> &pw_contact_pair_candidates,
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      &pw_pairs_in_contact);
};

#endif /* PWFINESEARCH_H_ */
