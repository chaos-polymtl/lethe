// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/contact_info.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_fine_search.h>

#include <deal.II/particles/particle.h>

#include <boost/range/adaptor/map.hpp>

#include <unordered_map>

using namespace dealii;

template <int dim, ContactType contact_type>
void
particle_particle_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
                                    &contact_pair_candidates,
  const double                       neighborhood_threshold,
  const std::vector<Tensor<1, dim>> &combined_offsets)
{
  // First iterating over adjacent_particles
  for (auto &&adjacent_particles_list :
       adjacent_particles | boost::adaptors::map_values)
    {
      if (adjacent_particles_list.empty())
        continue;

      // Gather information about particle 1
      auto &particle_one = adjacent_particles_list.begin()->second.particle_one;
      Point<dim, double> particle_one_location = particle_one->get_location();

      if constexpr (contact_type == ContactType::local_particle_particle ||
                    contact_type == ContactType::ghost_particle_particle)
        {
          // Iterating over each map which contains the contact information
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list.begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list.end();)
            {
              // Getting contact information and particle 2 as local
              // variables
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
                    adjacent_particles_list.erase(
                      adjacent_particles_list_iterator);
                }
              else
                {
                  // Save a translation that falls within the threshold
                  // Tensor<1, 3> offset_3d;
                  // for (int d = 0; d < dim; ++d)
                  //  offset_3d[d] = nearest_translation[d];

                  // adjacent_pair_information.periodic_offset = offset_3d;

                  ++adjacent_particles_list_iterator;
                }
            }
        }

      if constexpr (contact_type ==
                      ContactType::local_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_local_periodic_particle_particle)
        {
          // Iterating over each map which contains the contact information
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list.begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list.end();)
            {
              // Getting contact information and particle 2 as local
              // variables
              auto &adjacent_pair_information =
                adjacent_particles_list_iterator->second;
              auto &particle_two = adjacent_pair_information.particle_two;

              // Finding the properties of the particles in contact
              Point<dim, double> particle_two_real_location =
                particle_two->get_location();

              // Find minimum periodic distance squared between particle one and
              // the set of periodic images of particle two. The nearest
              // periodic image corresponds to a given periodic offset,
              // indicating if the particles are in contact through a periodic
              // edge, face, or corner.
              double min_square_distance = std::numeric_limits<double>::max();
              Tensor<1, dim> nearest_translation;

              for (const auto &translation : combined_offsets)
                {
                  // Check particle 1 against every periodic image of particle 2
                  double current_square_distance =
                    particle_one_location.distance_square(
                      particle_two_real_location + translation);

                  if (current_square_distance < min_square_distance)
                    {
                      min_square_distance = current_square_distance;
                      nearest_translation = translation;
                    }
                }

              // If simulation is well defined, there should be at most one
              // translation from combined_offsets that brings a periodic image
              // within a neighborhood threshold
              if (min_square_distance > neighborhood_threshold)
                {
                  adjacent_particles_list_iterator =
                    adjacent_particles_list.erase(
                      adjacent_particles_list_iterator);
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
          Point<dim, double> particle_two_real_location =
            particle_two->get_location();

          // Finding distance
          double min_square_distance = std::numeric_limits<double>::max();
          Tensor<1, dim> nearest_translation;

          for (const auto &translation : combined_offsets)
            {
              double current_square_distance =
                particle_one_location.distance_square(
                  particle_two_real_location + translation);

              if (current_square_distance < min_square_distance)
                {
                  min_square_distance = current_square_distance;
                  nearest_translation = translation;
                }
            }

          // If particles are within the threshold
          if (min_square_distance < neighborhood_threshold)
            {
              // Save a translation that falls within the threshold
              Tensor<1, 3> offset_3d;
              for (int d = 0; d < dim; ++d)
                offset_3d[d] = nearest_translation[d];

              auto &particle_one_contact_list =
                adjacent_particles[particle_one_id];

              particle_one_contact_list.emplace(
                particle_two_id,
                particle_particle_contact_info<dim>{particle_one,
                                                    particle_two,
                                                    Tensor<1, 3>(),
                                                    Tensor<1, 3>(),
                                                    offset_3d});
            }
        }
    }
}

template void
particle_particle_fine_search<2, ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 2>> &combined_offsets = {Tensor<1, 2>()});

template void
particle_particle_fine_search<2, ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 2>> &combined_offsets = {Tensor<1, 2>()});

template void
particle_particle_fine_search<2, ContactType::local_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 2>> &combined_offsets = {Tensor<1, 2>()});

template void
particle_particle_fine_search<2, ContactType::ghost_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 2>> &combined_offsets = {Tensor<1, 2>()});

template void
particle_particle_fine_search<
  2,
  ContactType::ghost_local_periodic_particle_particle>(
  typename DEM::dem_data_structures<2>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 2>> &combined_offsets = {Tensor<1, 2>()});

template void
particle_particle_fine_search<3, ContactType::local_particle_particle>(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 3>> &combined_offsets = {Tensor<1, 3>()});

template void
particle_particle_fine_search<3, ContactType::ghost_particle_particle>(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 3>> &combined_offsets = {Tensor<1, 3>()});
template void
particle_particle_fine_search<3, ContactType::local_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 3>> &combined_offsets = {Tensor<1, 3>()});
template void
particle_particle_fine_search<3, ContactType::ghost_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 3>> &combined_offsets = {Tensor<1, 3>()});
template void

particle_particle_fine_search<3, ContactType::ghost_local_periodic_particle_particle>(
  typename DEM::dem_data_structures<3>::particle_index_iterator_map const
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
                                  &contact_pair_candidates,
  const double                     neighborhood_threshold,
  const std::vector<Tensor<1, 3>> &combined_offsets = {Tensor<1, 3>()});
