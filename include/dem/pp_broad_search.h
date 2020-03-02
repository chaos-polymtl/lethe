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
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/pp_contact_info_struct.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef PPBROADSEARCH_H_
#  define PPBROADSEARCH_H_

/**
 * This class is used for broad particle-particle contact search. Broad search
 * is used to obtain all the particle-particle pairs in adjacent cells.
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class PPBroadSearch
{
public:
  PPBroadSearch<dim>();

  /**
   * Finds a vector of pairs (particle pairs) which shows the candidate
   * particle-particle collision pairs. These collision pairs will be used in
   * the fine search to investigate if they are in contact or not.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param cellNeighborList This vector is the output of FindCellNeighbors
   * class and shows the neighbor cells of each cell in the triangulation
   * @return A vector of pairs which contains all the particle pairs in adjacent
   * cells which are collision candidates
   */

  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
  find_PP_Contact_Pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
      &cellNeighborList);
};

#endif /* PPBROADSEARCH_H_ */
