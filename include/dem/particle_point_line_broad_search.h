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

#include <dem/data_containers.h>
#include <dem/disable_contacts.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

using namespace dealii;

#ifndef particle_point_line_broad_search_h
#  define particle_point_line_broad_search_h

/**
 * This class is used for broad particle-line and particle-point contact
 * search.This broad search is used to obtain all the particles located at
 * boundary cells with lines and points
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticlePointLineBroadSearch
{
public:
  ParticlePointLineBroadSearch();

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

  typename DEM::dem_data_structures<dim>::particle_point_candidates
  find_particle_point_contact_pairs(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::unordered_map<
      std::string,
      std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
      &boundary_cells_with_points);

  typename DEM::dem_data_structures<dim>::particle_point_candidates
  find_particle_point_contact_pairs(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::unordered_map<
      std::string,
      std::pair<typename Triangulation<dim>::active_cell_iterator, Point<dim>>>
                               &boundary_cells_with_points,
    const DisableContacts<dim> &disable_contacts_object);

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

  typename DEM::dem_data_structures<dim>::particle_line_candidates
  find_particle_line_contact_pairs(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::unordered_map<
      std::string,
      std::tuple<typename Triangulation<dim>::active_cell_iterator,
                 Point<dim>,
                 Point<dim>>> &boundary_cells_with_lines);

  typename DEM::dem_data_structures<dim>::particle_line_candidates
  find_particle_line_contact_pairs(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::unordered_map<
      std::string,
      std::tuple<typename Triangulation<dim>::active_cell_iterator,
                 Point<dim>,
                 Point<dim>>>  &boundary_cells_with_lines,
    const DisableContacts<dim> &disable_contacts_object);
};

#endif /* particle_point_line_broad_search_h */
