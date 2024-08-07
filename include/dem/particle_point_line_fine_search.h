/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_particle_point_line_fine_search_h
#define lethe_particle_point_line_fine_search_h

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/particle_point_line_contact_info.h>

#include <deal.II/particles/particle_handler.h>

#include <iostream>
#include <vector>

using namespace dealii;

/**
 * @brief Iterate over a map of pairs (particle_point_contact_candidates)
 * to see if the particle-point pairs are in contact or not. If they are in
 * contact, the normal overlap, normal vector of contact and contact normal
 * relative velocity are stored in a map which is the output of this function.
 *
 * @param particle_point_contact_candidates The output of particle-point broad
 * search which shows contact pair candidates.
 * @param neighborhood_threshold A value which defines the neighbor particles.
 * @param particle_point_pairs_in_contact A map which contains all the required
 * information (normal overlap, normal vector and contact normal relative
 * velocity) for calculation of particle-point contact force.
 */
template <int dim>
void
particle_point_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<dim>::particle_point_in_contact
    &particle_point_pairs_in_contact);

/**
 * @brief Iterate over a map of tuples (particle_line_contact_candidates) to
 * see if the particle-line pairs are in contact or not. If they are in
 * contact, the normal overlap, normal vector of contact and contact normal
 * relative velocity are stored in a map which is the output of this function.
 *
 * @param particle_line_contact_candidates The output of particle-line broad
 * search which shows contact pair candidates.
 * @param neighborhood_threshold A value which defines the neighbor particles.
 * @param particle_line_pairs_in_contact A map which contains all the required
 * information (normal overlap, normal vector and contact normal relative
 * velocity) for calculation of particle-line contact force.
 */
template <int dim>
void
particle_line_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<dim>::particle_line_in_contact
    &particle_line_pairs_in_contact);

/**
 * @brief This private function is used to find the projection of point_p on
 * a line with beginning and ending vertices of point_a and point_b,
 * respectively
 *
 * @param point_p A point which is going to be projected on the line
 * @param point_a Beginning point of the line
 * @param point_b Ending point of the line
 *
 * @return The projection of point_p on the line (from point_a to point_b)
 */
inline Point<3>
find_projection_point(const Point<3> &point_p,
                      const Point<3> &point_a,
                      const Point<3> &point_b);

#endif
