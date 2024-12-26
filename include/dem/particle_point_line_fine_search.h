// SPDX-FileCopyrightText: Copyright (c) 2020, 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_point_line_fine_search_h
#define lethe_particle_point_line_fine_search_h

#include <dem/data_containers.h>

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
template <int dim, DEM::SolverType solver_type>
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
template <int dim, DEM::SolverType solver_type>
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
