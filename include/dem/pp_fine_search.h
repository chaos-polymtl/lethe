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

#include "dem/dem_properties.h"
#include "dem/pp_contact_info_struct.h"
#include <deal.II/base/tensor.h>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <iostream>
#include <vector>

using namespace dealii;

#ifndef PPFINESEARCH_H_
#define PPFINESEARCH_H_

/**
 * This class is used for fine particle-particle contact search. Fine search
 * is used to find all the particle pairs which are physically in contact and
 * obtain all the required information for calculation of the contact force
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim, int spacedim = dim> class PPFineSearch {
public:
  PPFineSearch<dim, spacedim>();

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
   * @param contact_pair_candidates The output of broad search which shows
   * contact pair candidates
   * @param pairs_in_contact A vector of maps which stores all the particle
   * pairs which are physically in contact
   * @param pairs_in_contact_info A vector of maps which stores all the required
   * information for calculation of the contact force
   * @param dt DEM time step
   * @param particle_handler The particle handler of particles in the fine
   * search
   */

  void pp_Fine_Search(
      std::vector<std::pair<Particles::ParticleIterator<dim, spacedim>,
                            Particles::ParticleIterator<dim, spacedim>>>
          &contact_pair_candidates,
      std::vector<std::map<int, Particles::ParticleIterator<dim, spacedim>>>
          &pairs_in_contact,
      std::vector<std::map<int, pp_contact_info_struct<dim, spacedim>>>
          &pairs_in_contact_info,
      double dt, Particles::ParticleHandler<dim, spacedim> &particle_handler);
};

#endif /* PPFINESEARCH_H_ */
