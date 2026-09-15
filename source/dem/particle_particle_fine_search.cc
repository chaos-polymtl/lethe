// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_fine_search.h>

#include <deal.II/particles/particle.h>

#include <boost/range/adaptor/map.hpp>

#include <cmath>
#include <limits>
#include <unordered_map>

using namespace dealii;

template <int dim, ContactType contact_type>
void
particle_particle_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                                  &particle_container,
  adjacent_pairs_for_contact_t<dim, contact_type> &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                       &contact_pair_candidates,
  const double          neighborhood_threshold,
  const Tensor<1, dim> &periodic_offset_per_direction)
{
  // First iterating over adjacent_particles
  for (auto &&adjacent_particles_list :
       adjacent_particles | boost::adaptors::map_values)
    {
      auto &second_particles = adjacent_particles_list.second_particles;

      if (second_particles.empty())
        continue;

      // Gather information about particle 1. The iterator to particle 1 is
      // shared by all its contacts and stored once in the adjacency value.
      auto              &particle_one = adjacent_particles_list.particle_one;
      Point<dim, double> particle_one_location = particle_one->get_location();

      // For non-periodic contacts
      if constexpr (contact_type == local_particle_particle ||
                    contact_type == ghost_particle_particle)
        {
          // Iterating over each map which contains the contact information
          for (auto adjacent_particles_list_iterator = second_particles.begin();
               adjacent_particles_list_iterator != second_particles.end();)
            {
              // Getting contact information and particle 2 as local variables
              auto &adjacent_pair_information =
                adjacent_particles_list_iterator->second;
              auto &particle_two = adjacent_pair_information.particle_two;

              // Finding the properties of the particles in contact
              Point<dim, double> particle_two_location =
                particle_two->get_location();

              const double square_distance =
                particle_one_location.distance_square(particle_two_location);
              if (square_distance > neighborhood_threshold)
                {
                  adjacent_particles_list_iterator =
                    second_particles.erase(adjacent_particles_list_iterator);
                }
              else
                {
                  ++adjacent_particles_list_iterator;
                }
            }
        }

      // For periodic contacts
      if constexpr (contact_type == local_periodic_particle_particle ||
                    contact_type == ghost_periodic_particle_particle ||
                    contact_type == ghost_local_periodic_particle_particle)
        {
          // Iterating over each map which contains the contact information
          for (auto adjacent_particles_list_iterator = second_particles.begin();
               adjacent_particles_list_iterator != second_particles.end();)
            {
              // Getting contact information and particle 2 as local variables
              auto &adjacent_pair_information =
                adjacent_particles_list_iterator->second;
              auto &particle_two = adjacent_pair_information.particle_two;

              // Finding the properties of the particles in contact
              Point<dim, double> particle_two_real_location =
                particle_two->get_location();

              // Find the translation that brings the nearest periodic image of
              // particle two to particle one. Since periodic directions are
              // axis-aligned and independent, this translation is found
              // directly, one direction at a time (minimum image convention),
              // instead of searching over every combination of periodic
              // offsets. A pair whose nearest image requires no translation on
              // any periodic direction is not a periodic contact (it is
              // already handled by the non-periodic contact types), so it is
              // treated as out of range.
              Tensor<1, dim> nearest_translation;
              bool           found_periodic_translation = false;

              for (int d = 0; d < dim; ++d)
                {
                  if (periodic_offset_per_direction[d] != 0.)
                    {
                      const double delta = particle_one_location[d] -
                                           particle_two_real_location[d];
                      nearest_translation[d] =
                        std::round(delta / periodic_offset_per_direction[d]) *
                        periodic_offset_per_direction[d];
                      found_periodic_translation |=
                        (nearest_translation[d] != 0.);
                    }
                }

              const double min_square_distance =
                found_periodic_translation ?
                  particle_one_location.distance_square(
                    particle_two_real_location + nearest_translation) :
                  std::numeric_limits<double>::max();

              // If simulation is well defined, there should be at most one
              // periodic image that brings the particles within a
              // neighborhood threshold.
              if (min_square_distance > neighborhood_threshold)
                {
                  adjacent_particles_list_iterator =
                    second_particles.erase(adjacent_particles_list_iterator);
                }
              else
                {
                  // Save a translation that falls within the threshold
                  Tensor<1, 3> offset_3d;
                  for (int d = 0; d < dim; ++d)
                    offset_3d[d] = nearest_translation[d];

                  adjacent_pair_information.periodic_offset = offset_3d;

                  ++adjacent_particles_list_iterator;
                }
            }
        }
    }

  // Now iterating over contact_pair_candidates (maps of pairs), which
  // is the output of broad search. If a pair is in vicinity (distance <
  // threshold), it is added to the adjacent_particles
  for (auto &[particle_one_id, second_particle_container] :
       contact_pair_candidates)
    {
      if (second_particle_container.empty())
        continue;

      auto               particle_one = particle_container.at(particle_one_id);
      Point<dim, double> particle_one_location = particle_one->get_location();

      for (const types::particle_index &particle_two_id :
           second_particle_container)
        {
          auto particle_two = particle_container.at(particle_two_id);

          // For non-periodic contacts
          if constexpr (contact_type == local_particle_particle ||
                        contact_type == ghost_particle_particle)
            {
              Point<dim, double> particle_two_location =
                particle_two->get_location();

              const double square_distance =
                particle_one_location.distance_square(particle_two_location);

              if (square_distance < neighborhood_threshold)
                {
                  auto &particle_one_contact_list =
                    adjacent_particles[particle_one_id];
                  particle_one_contact_list.particle_one = particle_one;

                  particle_one_contact_list.second_particles.emplace(
                    particle_two_id,
                    particle_particle_contact_info<dim>{particle_two,
                                                        Tensor<1, 3>(),
                                                        Tensor<1, 3>()});
                }
            }

          // For periodic contacts
          if constexpr (contact_type == local_periodic_particle_particle ||
                        contact_type == ghost_periodic_particle_particle ||
                        contact_type == ghost_local_periodic_particle_particle)
            {
              Point<dim, double> particle_two_real_location =
                particle_two->get_location();

              // Find the translation that brings the nearest periodic image of
              // particle two to particle one. Since periodic directions are
              // axis-aligned and independent, this translation is found
              // directly, one direction at a time (minimum image convention),
              // instead of searching over every combination of periodic
              // offsets. A pair whose nearest image requires no translation on
              // any periodic direction is not a periodic contact (it is
              // already handled by the non-periodic contact types), so it is
              // treated as out of range.
              Tensor<1, dim> nearest_translation;
              bool           found_periodic_translation = false;

              for (int d = 0; d < dim; ++d)
                {
                  if (periodic_offset_per_direction[d] != 0.)
                    {
                      const double delta = particle_one_location[d] -
                                           particle_two_real_location[d];
                      nearest_translation[d] =
                        std::round(delta / periodic_offset_per_direction[d]) *
                        periodic_offset_per_direction[d];
                      found_periodic_translation |=
                        (nearest_translation[d] != 0.);
                    }
                }

              const double min_square_distance =
                found_periodic_translation ?
                  particle_one_location.distance_square(
                    particle_two_real_location + nearest_translation) :
                  std::numeric_limits<double>::max();

              // If simulation is well defined, there should be at most one
              // periodic image that brings the particles within a
              // neighborhood threshold.
              if (min_square_distance < neighborhood_threshold)
                {
                  // Save a translation that falls within the threshold
                  Tensor<1, 3> offset_3d;
                  for (int d = 0; d < dim; ++d)
                    offset_3d[d] = nearest_translation[d];

                  auto &particle_one_contact_list =
                    adjacent_particles[particle_one_id];
                  particle_one_contact_list.particle_one = particle_one;

                  particle_one_contact_list.second_particles.emplace(
                    particle_two_id,
                    periodic_particle_particle_contact_info<dim>{
                      {particle_two, Tensor<1, 3>(), Tensor<1, 3>()},
                      offset_3d});
                }
            }
        }
    }
}

// 2D templates
template void
particle_particle_fine_search<2, local_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 2> &periodic_offset_per_direction);


// 3D templates
template void
particle_particle_fine_search<3, local_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                     &contact_pair_candidates,
  const double        neighborhood_threshold,
  const Tensor<1, 3> &periodic_offset_per_direction);
