// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/contact_info.h>
#include <dem/particle_point_line_fine_search.h>

#include <deal.II/particles/particle_handler.h>
using namespace dealii;

// In this function, the output of particle-point broad search is investigated
// to see if the pairs are in contact or not. If they are in contact, the normal
// overlap, normal vector and normal relative velocity of the contact are
// calculated. The output of this function is used for calculation of the
// contact force
template <int dim, typename PropertiesIndex>
void
particle_point_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<dim>::particle_point_in_contact
    &particle_point_pairs_in_contact)
{
  // Iterating over contact candidates from broad search. If a particle-point
  // pair is in contact (distance > 0) it is inserted into the output of this
  // function (particle_point_pairs_in_contact)
  for (auto contact_pair_candidates_iterator =
         particle_point_contact_candidates.begin();
       contact_pair_candidates_iterator !=
       particle_point_contact_candidates.end();
       ++contact_pair_candidates_iterator)
    {
      // Get the value of the map (pair candidate) from the
      // contact_pair_candidates_iterator
      const particle_point_contact_info<dim> &pair_candidates =
        contact_pair_candidates_iterator->second;

      // Get the particle, particle diameter and boundary vertex location once
      Particles::ParticleIterator<dim> particle = pair_candidates.particle;
      double                           particle_diameter =
        particle->get_properties()[PropertiesIndex::dp];

      Point<3> vertex_location   = pair_candidates.point;
      Point<3> particle_location = [&] {
        if constexpr (dim == 3)
          return particle->get_location();
        else
          return point_nd_to_3d(particle->get_location());
      }();

      // Calculation of the square_distance between the particle and boundary
      // vertex
      const double square_distance =
        (particle_diameter / 2.0) -
        vertex_location.distance_square(particle_location);

      // If the distance is larger than neighborhood threshold, then the
      // particle-point pair are in contact
      if (square_distance > neighborhood_threshold)
        {
          // Adding contact info to the sample
          particle_point_pairs_in_contact.emplace(particle->get_id(),
                                                  pair_candidates);
        }
    }
}

// In this function, the output of particle-line broad search is investigated
// to see if the pairs are in contact or not. If they are in contact, the normal
// overlap, normal vector and normal relative velocity of the contact are
// calculated. The output of this function is used for calculation of the
// contact force
template <int dim, typename PropertiesIndex>
void
particle_line_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<dim>::particle_line_in_contact
    &particle_line_pairs_in_contact)
{
  // Iterating over contact candidates from broad search. If a particle-line
  // pair is in contact (distance > 0) it is inserted into the output of this
  // function (particle_line_pairs_in_contact);
  for (auto contact_pair_candidates_iterator =
         particle_line_contact_candidates.begin();
       contact_pair_candidates_iterator !=
       particle_line_contact_candidates.end();
       ++contact_pair_candidates_iterator)
    {
      // Get the value of the map (pair candidate) from the
      // contact_pair_candidates_iterator
      const particle_line_contact_info<dim> &pair_candidates =
        contact_pair_candidates_iterator->second;

      // Get the particle, particle diameter and the locations of beginning
      // and ending vertices of the boundary line
      Particles::ParticleIterator<dim> particle = pair_candidates.particle;
      double                           particle_diameter =
        particle->get_properties()[PropertiesIndex::dp];

      Point<3> vertex_one_location = pair_candidates.point_one;
      Point<3> vertex_two_location = pair_candidates.point_two;
      Point<3> particle_location   = [&] {
        if constexpr (dim == 3)
          return particle->get_location();
        else
          return point_nd_to_3d(particle->get_location());
      }();

      // For finding the particle-line distance, the projection of the particle
      // on the line should be obtained
      Point<3> projection = find_projection_point(particle_location,
                                                  vertex_one_location,
                                                  vertex_two_location);

      // Calculation of the distance between the particle and boundary line
      const double square_distance =
        (particle_diameter / 2.0) -
        projection.distance_square(particle_location);

      // If the distance is positive, then the particle-line pair are in
      // contact
      if (square_distance > neighborhood_threshold)
        {
          // Creating a sample from the particle_point_line_contact_info
          // and adding contact info to the sample
          particle_line_pairs_in_contact.emplace(particle->get_id(),
                                                 pair_candidates);
        }
    }
}

Point<3>
find_projection_point(const Point<3> &point_p,
                      const Point<3> &point_a,
                      const Point<3> &point_b)
{
  Tensor<1, 3> vector_ab = point_b - point_a;
  Tensor<1, 3> vector_ap = point_p - point_a;

  Point<3> projection =
    point_a + ((vector_ap * vector_ab) / (vector_ab * vector_ab)) * vector_ab;

  return projection;
}

template void
particle_point_fine_search<2, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_point_fine_search<3, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_line_fine_search<2, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_line_in_contact
    &particle_line_pairs_in_contact);

template void
particle_line_fine_search<3, DEM::DEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_line_in_contact
    &particle_line_pairs_in_contact);

template void
particle_point_fine_search<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_point_fine_search<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_line_fine_search<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_line_in_contact
    &particle_line_pairs_in_contact);

template void
particle_line_fine_search<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_line_in_contact
    &particle_line_pairs_in_contact);

template void
particle_point_fine_search<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_point_fine_search<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_point_candidates
              &particle_point_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_point_in_contact
    &particle_point_pairs_in_contact);

template void
particle_line_fine_search<2, DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<2>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<2>::particle_line_in_contact
    &particle_line_pairs_in_contact);

template void
particle_line_fine_search<3, DEM::DEMMPProperties::PropertiesIndex>(
  const DEM::dem_data_structures<3>::particle_line_candidates
              &particle_line_contact_candidates,
  const double neighborhood_threshold,
  DEM::dem_data_structures<3>::particle_line_in_contact
    &particle_line_pairs_in_contact);
