// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_fine_search.h>

#include <deal.II/particles/particle.h>

#include <boost/range/adaptor/map.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <utility>

using namespace dealii;

template <int dim, ContactType contact_type>
void
particle_particle_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                                  &particle_container,
  adjacent_pairs_for_contact_t<dim, contact_type> &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                                &contact_pair_candidates,
  const double                   neighborhood_threshold,
  const std::array<double, dim> &periodic_offset_per_direction)
{
  // First iterating over adjacent_particles
  for (auto &&adjacent_particles_list :
       adjacent_particles | boost::adaptors::map_values)
    {
      // List of potential particles
      auto &second_particles = adjacent_particles_list.second_particles;

      // If the list is empty, we continue.
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
              // Getting the contact information and particle 2 as local
              // variables
              auto &adjacent_pair_information =
                adjacent_particles_list_iterator->second;
              auto &particle_two = adjacent_pair_information.particle_two;

              // Finding the properties of particles 2.
              Point<dim, double> particle_two_real_location =
                particle_two->get_location();

              // Reuse the periodic offset found on a previous call: for a
              // persisting contact this translation practically never
              // changes between fine search calls. This avoids the
              // minimum periodic image search with the signed periodic
              // distance.
              Tensor<1, dim> cached_translation =
                adjacent_pair_information.periodic_offset;
              const double cached_square_distance =
                particle_one_location.distance_square(
                  particle_two_real_location + cached_translation);

              // Using the caches periodic offset, if the shared distance
              // respect the neighborhood threshold, this means that the
              // particles stayed on their respective side of the PBC.
              if (cached_square_distance <= neighborhood_threshold)
                ++adjacent_particles_list_iterator;

              // If not, either one particle crossed the PBC or both crossed it
              // or they just moved apart.
              else
                {
                  // Cached periodic offset no longer holds. Fall back to the
                  // distance computation using the minimum image convention.
                  // If the distance between the particles is higher than the
                  // neighbor threshold, we erase the potential contact from the
                  // list.
                  // The nearest_translation is the periodic offset.
                  const auto [nearest_translation, found_periodic_translation] =
                    nearest_periodic_translation<dim>(
                      particle_one_location,
                      particle_two_real_location,
                      periodic_offset_per_direction);

                  // If a periodic offset was found, the min_square_distance is
                  // computed. Otherwise, we know that the particle are far
                  // apart, thus we impose the min_square_distance to a high
                  // value (max())
                  const double min_square_distance =
                    found_periodic_translation ?
                      particle_one_location.distance_square(
                        particle_two_real_location + nearest_translation) :
                      std::numeric_limits<double>::max();

                  // If the min_square_distance does not respect the
                  // neighborhood threshold at this point, this means that the
                  // particle officially moved apart. We can remove the
                  // potential contact from the list.
                  if (min_square_distance > neighborhood_threshold)
                    adjacent_particles_list_iterator =
                      second_particles.erase(adjacent_particles_list_iterator);

                  // Otherwise, this means that one or both particle crossed the
                  // PBC. We update the periodic_offset in the periodic
                  // contact_info.
                  else
                    {
                      adjacent_pair_information.periodic_offset =
                        nearest_translation;
                      ++adjacent_particles_list_iterator;
                    }
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

              // nearest_periodic_translation returns periodic translation that
              // brings two particle the closest together. (if they are close
              // enough)
              const auto [nearest_translation, found_periodic_translation] =
                nearest_periodic_translation<dim>(
                  particle_one_location,
                  particle_two_real_location,
                  periodic_offset_per_direction);

              // If is possible that the two particle moved far appart so much
              // that the round operation in the nearest_periodic_translation
              // function return a 0. In this case, found_periodic_translation
              // will be at false. Thus, we need to make sure that this isn't
              // the case by imposing min_square_distance to max().
              const double min_square_distance =
                found_periodic_translation ?
                  particle_one_location.distance_square(
                    particle_two_real_location + nearest_translation) :
                  std::numeric_limits<double>::max();

              // If the neighborhood_threshold is respected, we add particle 2
              // to particle 1 potential contact list. Otherwise, we do nothing.
              if (min_square_distance < neighborhood_threshold)
                {
                  auto &particle_one_contact_list =
                    adjacent_particles[particle_one_id];
                  particle_one_contact_list.particle_one = particle_one;

                  particle_one_contact_list.second_particles.emplace(
                    particle_two_id,
                    periodic_particle_particle_contact_info<dim>{
                      {particle_two, Tensor<1, 3>(), Tensor<1, 3>()},
                      nearest_translation});
                }
            }
        }
    }
}

template <int dim>
std::pair<Tensor<1, dim>, bool>
nearest_periodic_translation(
  const Point<dim, double>      &particle_one_location,
  const Point<dim, double>      &particle_two_real_location,
  const std::array<double, dim> &periodic_offset_per_direction)
{
  Tensor<1, dim> nearest_translation;
  bool           found_periodic_translation = false;

  // Loop on every direction
  for (int d = 0; d < dim; ++d)
    {
      // If the periodic offset associated with this direction is zero,
      // this means that the direction is not periodic, thus we skip it.
      // Also, it prevents a division by zero a bit later in the code.
      if (periodic_offset_per_direction[d] != 0.)
        {
          // Real distance between p1 and p2
          const double delta =
            particle_one_location[d] - particle_two_real_location[d];

          // If particle 1 (p1) is on one side of the triangulation and p2 is
          // on the other side and both sides/boundaries are linked through a
          // PBC, the distance between those two particle should
          // be around one period (~ 0.9), thus when we round, we get 1.0,
          // which we then multiply to the actual periodic_offset. When p1 and
          // p2 and near periodic corner, this method will also work because the
          // rounding will give 1.0 for more than one direction.
          nearest_translation[d] =
            std::round(delta / periodic_offset_per_direction[d]) *
            periodic_offset_per_direction[d];

          // If the nearest_translation got rounded to 0., this means that
          // the particles are in the same cell or in a real neighboring cell
          // or that the current direction is not linking this periodic
          // contact. (By real neighboring cell, we mean a neighboring cell that
          // is not periodic neighboring cell)
          found_periodic_translation |= (nearest_translation[d] != 0.);
        }
    }
  return {nearest_translation, found_periodic_translation};
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
  const double                 neighborhood_threshold,
  const std::array<double, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 2> &periodic_offset_per_direction);

template void
particle_particle_fine_search<2, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 2> &periodic_offset_per_direction);


// 3D templates
template void
particle_particle_fine_search<3, local_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 3> &periodic_offset_per_direction);

template void
particle_particle_fine_search<3, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                              &contact_pair_candidates,
  const double                 neighborhood_threshold,
  const std::array<double, 3> &periodic_offset_per_direction);
