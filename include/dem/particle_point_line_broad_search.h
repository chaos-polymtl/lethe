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

#include <deal.II/distributed/tria.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef PARTICLEPOINTLINEBROADSEARCH_H_
#define PARTICLEPOINTLINEBROADSEARCH_H_

/**
 * This class is used for broad particle-line and particle-point contact
 * search.This broad search is used to obtain all the particles located at
 * boundary cells with lines and points
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim> class ParticlePointLineBroadSearch {
public:
  ParticlePointLineBroadSearch<dim>();

  /**
   * Finds a map of pairs (pair of particle and the boundary vertex location)
   * which shows the candidate particle-point collision pairs. These collision
   * pairs will be investigated in the fine search to check if they are in
   * contact or not
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param boundary_cells_with_points A container of cells which are located at
   * boundaries with only one vertex
   * @return A map of pairs. Each element of map (pair) contains a contact pair
   * (particle located near boundaries with vertices and the vertex location)
   */

  std::map<int, std::pair<Particles::ParticleIterator<dim>, Point<dim>>>
  find_Particle_Point_Contact_Pairs(
      const Particles::ParticleHandler<dim> &particle_handler,
      const std::vector<std::pair<
          typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
          &boundary_cells_with_points);

  /**
   * Finds a map of tuples (tuple of particle and the locations of beginning and
   * ending vertices of the boundary lines) which shows the candidate
   * particle-line collision pairs. These collision pairs will be investigated
   * in the fine search to check if they are in contact or not
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param boundary_cells_with_lines A container of cells which are located at
   * boundaries with only one line
   * @return A map of tuples. Each element of map (tuple) contains a particle
   * and the locations of beginning and ending vertices of the boundary lines
   */

  std::map<int,
           std::tuple<Particles::ParticleIterator<dim>, Point<dim>, Point<dim>>>
  find_Particle_Line_Contact_Pairs(
      const Particles::ParticleHandler<dim> &particle_handler,
      const std::vector<
          std::tuple<typename Triangulation<dim>::active_cell_iterator,
                     Point<dim>, Point<dim>>> &boundary_cells_with_lines);
};

#endif /* PARTICLEPOINTLINEBROADSEARCH_H_ */
