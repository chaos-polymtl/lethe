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
   * @param cells_local_neighbor_list This vector is the output of
   * find_cell_neighbors class and shows the local neighbor cells of all local
   * celsl in the triangulation
   * @param cells_ghost_neighbor_list This vector is the output of
   * find_cell_neighbors class and shows the ghost neighbor cells of all local
   * celsl in the triangulation
   * @param local_contact_pair_candidates A map of vectors which contains all
   * the local-local particle pairs in adjacent cells which are collision
   * candidates
   * @param ghost_contact_pair_candidates A map of vectors which contains all
   * the local-ghost particle pairs in adjacent cells which are collision
   * candidates
   */

  void
  find_PP_Contact_Pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      *cells_local_neighbor_list,
    const std::vector<
      std::vector<typename Triangulation<dim>::active_cell_iterator>>
      *                              cells_ghost_neighbor_list,
    std::map<int, std::vector<int>> &local_contact_pair_candidates,
    std::map<int, std::vector<int>> &ghost_contact_pair_candidates);
};

#endif /* PPBROADSEARCH_H_ */
