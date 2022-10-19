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
#include <dem/particle_point_line_contact_info_struct.h>

#include <deal.II/particles/particle_handler.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_point_line_fine_search_h
#  define particle_point_line_fine_search_h

/**
 * This class is used for fine particle-point and particle-line contact search.
 * Fine search is used to find all the particles which are physically in contact
 * with boundary lines and points, and obtain all the required information for
 * calculation of the corresponding contact force
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticlePointLineFineSearch
{
public:
  ParticlePointLineFineSearch<dim>();

  /**
   * Iterates over a map of pairs (particle_point_contact_candidates) to see if
   * the particle-point pairs are in contact or not. If they are in contact, the
   * normal overlap, normal vector of contact and contact normal relative
   * velocity are stored in a map which is the output of this function
   *
   * @param particle_point_contact_candidates The output of particle-point broad
   * search which shows contact pair candidates
   * @param neighborhood_threshold A value which defines the neighbor particles
   * @return A map which contains all the required information (normal overlap,
   * normal vector and contact normal relative velocity) for calculation of
   * particle-point contact force
   */

  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
  particle_point_fine_search(
    const typename DEM::dem_data_structures<dim>::particle_point_candidates
      &          particle_point_contact_candidates,
    const double neighborhood_threshold);

  /**
   * Iterates over a map of tuples (particle_line_contact_candidates) to see if
   * the particle-line pairs are in contact or not. If they are in contact, the
   * normal overlap, normal vector of contact and contact normal relative
   * velocity are stored in a map which is the output of this function
   *
   * @param particle_line_contact_candidates The output of particle-line broad
   * search which shows contact pair candidates
   * @param neighborhood_threshold A value which defines the neighbor particles
   * @return A map which contains all the required information (normal overlap,
   * normal vector and contact normal relative velocity) for calculation of
   * particle-line contact force
   */

  typename DEM::dem_data_structures<dim>::particle_point_line_contact_info
  particle_line_fine_search(
    const typename DEM::dem_data_structures<dim>::particle_line_candidates
      &          particle_line_contact_candidates,
    const double neighborhood_threshold);

private:
  /** This private function is used to find the projection of point_p on
   * a line with beginning and ending vertices of point_a and point_b,
   * respectively
   * @param point_p A point which is going to be projected on the line
   * @param point_a Beginning point of the line
   * @param point_b Ending point of the line
   * @return The projection of point_p on the line (from point_a to point_b)
   */

  Point<3>
  find_projection_point(const Point<3> &point_p,
                        const Point<3> &point_a,
                        const Point<3> &point_b);
};

#endif /* particle_point_line_fine_search_h */
