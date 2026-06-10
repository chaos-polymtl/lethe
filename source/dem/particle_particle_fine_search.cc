// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>

#include <dem/contact_info.h>
#include <dem/contact_type.h>
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
  adjacent_pairs_for_contact_t<dim, contact_type> &adjacent_particles,
  const typename DEM::dem_data_structures<dim>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<dim>::particle_index_tensor_map
                       &particle_index_periodic_offset_map,
  const Tensor<1, dim> &combined_periodic_offset)
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

      // For non-periodic contacts
      if constexpr (contact_type == local_particle_particle ||
                    contact_type == ghost_particle_particle)
        {
          // Iterating over each map that contains the contact information
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list.begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list.end();)
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
                    adjacent_particles_list.erase(
                      adjacent_particles_list_iterator);
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
          // Iterating over each map that contains the contact information
          for (auto adjacent_particles_list_iterator =
                 adjacent_particles_list.begin();
               adjacent_particles_list_iterator !=
               adjacent_particles_list.end();)
            {
              // Getting contact information and particle 2 as local variables
              auto &adjacent_pair_information =
                adjacent_particles_list_iterator->second;
              auto &particle_two_information =
                adjacent_pair_information.particle_two;

              // Finding the properties of the particles in contact
              Point<dim> particle_two_real_location =
                particle_two_information->get_location();

              Tensor<1, dim> periodic_offset_applied_to_particle_one =
                particle_index_periodic_offset_map[particle_one->get_id()];

              Tensor<1, dim> periodic_offset_applied_to_particle_two =
                particle_index_periodic_offset_map[particle_two_information
                                                     ->get_id()];

              // (Particle one location - Particle two location)
              double square_distance = particle_one_location.distance_square(
                particle_two_real_location +
                periodic_offset_applied_to_particle_one -
                combined_periodic_offset -
                periodic_offset_applied_to_particle_two);

              if (square_distance > neighborhood_threshold)
                {
                  adjacent_particles_list_iterator =
                    adjacent_particles_list.erase(
                      adjacent_particles_list_iterator);
                }
              else
                {
                  // Save a translation that falls within the threshold
                  if constexpr (dim == 2)
                    adjacent_pair_information.periodic_offset =
                      tensor_nd_to_3d(combined_periodic_offset);
                  else
                    {
                      adjacent_pair_information.periodic_offset =
                        combined_periodic_offset;
                    }
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
      // If this particle has no other particle close to it.
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

                  particle_one_contact_list.emplace(
                    particle_two_id,
                    particle_particle_contact_info<dim>{particle_one,
                                                        particle_two,
                                                        Tensor<1, 3>(),
                                                        Tensor<1, 3>()});
                }
            }

          // For periodic contacts
          if constexpr (contact_type == local_periodic_particle_particle ||
                        contact_type == ghost_periodic_particle_particle ||
                        contact_type == ghost_local_periodic_particle_particle)
            {
              Point<dim, double> particle_two_location =
                particle_two->get_location();

              Tensor<1, dim> periodic_offset_applied_to_particle_one =
                particle_index_periodic_offset_map[particle_one->get_id()];

              Tensor<1, dim> periodic_offset_applied_to_particle_two =
                particle_index_periodic_offset_map[particle_two->get_id()];

              std::cout << " Periodic offset: " << combined_periodic_offset << std::endl;
              Tensor<1, dim> total_offset =
                periodic_offset_applied_to_particle_one -
                combined_periodic_offset -
                periodic_offset_applied_to_particle_two;

              // (Particle one location - Particle two location)
              double square_distance = particle_one_location.distance_square(
                particle_two_location + total_offset);

              std::cout << " Particle 1-2 id    : " << particle_one->get_id()
                        << " : " << particle_two->get_id() << std::endl;
              std::cout << " Particle 1 location: " << particle_one_location
                        << std::endl;
              std::cout << " Particle 2 location: " << particle_two_location
                        << std::endl;
              std::cout << " Particle 1 applied periodic offset: "
                        << periodic_offset_applied_to_particle_one << std::endl;
              std::cout << " Particle 2 applied periodic offset: "
                        << periodic_offset_applied_to_particle_two << std::endl;
              std::cout << " Particle 2 relative position:" << particle_two_location + total_offset << std::endl;
              std::cout << " Periodic offset    : " << combined_periodic_offset
                        << std::endl;
              std::cout << " Square_distance    : " << square_distance
                        << std::endl;
              std::cout << " Threshold          : " << neighborhood_threshold
                        << std::endl;

              if (square_distance < neighborhood_threshold)
                {
                  std::cout << " True" << std::endl;
                  //  Save a translation that falls within the threshold
                  Tensor<1, 3> offset_3d = [total_offset]() -> Tensor<1, 3> {
                    if constexpr (dim == 2)
                      return tensor_nd_to_3d(total_offset);
                    else
                      return total_offset;
                  }();

                  auto &particle_one_contact_list =
                    adjacent_particles[particle_one_id];

                  particle_one_contact_list.emplace(
                    particle_two_id,
                    periodic_particle_particle_contact_info<dim>{
                      {particle_one,
                       particle_two,
                       Tensor<1, 3>(),
                       Tensor<1, 3>()},
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
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<2>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 2> &combined_periodic_offset);

template void
particle_particle_fine_search<2, ghost_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<2>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 2> &combined_periodic_offset);

template void
particle_particle_fine_search<2, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<2>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 2> &combined_periodic_offset);

template void
particle_particle_fine_search<2, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<2>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 2> &combined_periodic_offset);

template void
particle_particle_fine_search<2, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<2>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<2>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<2>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 2> &combined_periodic_offset);


// 3D templates
template void
particle_particle_fine_search<3, local_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<3>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 3> &combined_periodic_offset);

template void
particle_particle_fine_search<3, ghost_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<3>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 3> &combined_periodic_offset);

template void
particle_particle_fine_search<3, local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<3>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 3> &combined_periodic_offset);

template void
particle_particle_fine_search<3, ghost_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<3>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 3> &combined_periodic_offset);

template void
particle_particle_fine_search<3, ghost_local_periodic_particle_particle>(
  const typename DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container,
  typename DEM::dem_data_structures<3>::periodic_adjacent_particle_pairs
    &adjacent_particles,
  const typename DEM::dem_data_structures<3>::particle_particle_candidates
              &contact_pair_candidates,
  const double neighborhood_threshold,
  typename DEM::dem_data_structures<3>::particle_index_tensor_map
                     &particle_index_periodic_offset_map,
  const Tensor<1, 3> &combined_periodic_offset);
