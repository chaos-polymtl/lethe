// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/update_local_particle_containers.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

template <int dim>
void
update_particle_container(
  typename DEM::dem_data_structures<dim>::particle_index_iterator_map
                                        &particle_container,
  const Particles::ParticleHandler<dim> *particle_handler)
{
  // Clear the current particle container
  particle_container.clear();

  // Update the local particle container with local particles
  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      particle_container.emplace(particle_iterator->get_id(),
                                 particle_iterator);
    }

  // Update the local particle container with ghost particles
  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      particle_container.emplace(particle_iterator->get_id(),
                                 particle_iterator);
    }
}

template <int dim, typename pairs_structure, ContactType contact_type>
void
update_contact_container_iterators(
  pairs_structure &pairs_in_contact,
  const typename DEM::dem_data_structures<dim>::particle_index_iterator_map
    &particle_container)
{
  // Get trigger for clearing contact structures
  auto      *action_manager = DEMActionManager::get_action_manager();
  const bool clear_tangential_displacement =
    action_manager->check_clear_tangential_displacement();

  // Loop over particle-object (particle/wall/line/point) pairs in contact
  for (auto pairs_in_contact_iterator = pairs_in_contact.begin();
       pairs_in_contact_iterator != pairs_in_contact.end();)
    {
      // Get the adjacent objects content
      auto adjacent_pairs_content = &pairs_in_contact_iterator->second;

      if constexpr (contact_type == ContactType::local_particle_particle ||
                    contact_type == ContactType::ghost_particle_particle ||
                    contact_type ==
                      ContactType::local_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_periodic_particle_particle ||
                    contact_type ==
                      ContactType::ghost_local_periodic_particle_particle ||
                    contact_type == ContactType::particle_wall ||
                    contact_type == ContactType::particle_floating_wall)
        {
          // Get current particle id
          unsigned int particle_id = pairs_in_contact_iterator->first;

          // Find if particle one is still on this process
          auto particle_one_container = particle_container.find(particle_id);

          // Particle one is not on this process anymore, remove the container
          // and move on
          if (particle_one_container == particle_container.end())
            {
              pairs_in_contact_iterator =
                pairs_in_contact.erase(pairs_in_contact_iterator);
              continue;
            }

          // Loop over all the other objects of contact_type in contact
          for (auto adjacent_map_iterator = adjacent_pairs_content->begin();
               adjacent_map_iterator != adjacent_pairs_content->end();)
            {
              if constexpr (contact_type ==
                              ContactType::local_particle_particle ||
                            contact_type ==
                              ContactType::ghost_particle_particle ||
                            contact_type ==
                              ContactType::local_periodic_particle_particle ||
                            contact_type ==
                              ContactType::ghost_periodic_particle_particle ||
                            contact_type ==
                              ContactType::
                                ghost_local_periodic_particle_particle)
                {
                  unsigned int particle_two_id = adjacent_map_iterator->first;

                  // Find if particle two is still on this process
                  auto particle_two_container =
                    particle_container.find(particle_two_id);

                  // Particle two is not on this process anymore, remove the
                  // container and move on
                  if (particle_two_container == particle_container.end())
                    {
                      adjacent_map_iterator =
                        adjacent_pairs_content->erase(adjacent_map_iterator);
                      continue;
                    }

                  // Both particles are still on this process, update their
                  // iterators
                  adjacent_map_iterator->second.particle_one =
                    particle_one_container->second;
                  adjacent_map_iterator->second.particle_two =
                    particle_two_container->second;

                  if (clear_tangential_displacement)
                    {
                      adjacent_map_iterator->second.tangential_displacement
                        .clear();
                      adjacent_map_iterator->second
                        .rolling_resistance_spring_torque.clear();
                    }
                  ++adjacent_map_iterator;
                }

              if constexpr (contact_type == ContactType::particle_wall ||
                            contact_type == ContactType::particle_floating_wall)
                {
                  // Particle iterator is updated
                  adjacent_map_iterator->second.particle =
                    particle_container.at(particle_id);
                  ++adjacent_map_iterator;
                }
            }
        }

      if constexpr (contact_type == ContactType::particle_floating_mesh)
        {
          for (auto adjacent_map_iterator = adjacent_pairs_content->begin();
               adjacent_map_iterator != adjacent_pairs_content->end();
               ++adjacent_map_iterator)
            {
              // Get the particle in contact with the solid of the floating mesh
              unsigned int particle_id = adjacent_map_iterator->first;

              // Update the particle of the contact
              adjacent_map_iterator->second.particle =
                particle_container.at(particle_id);

              // Note : particle_triangle_cell_from_mesh_in_contact is stored as
              // <cell iterator, <particle id, particle-wall info>> and
              // particle_wall_in_contact as is stored as
              // <particle id, <face id, particle-wall info>>,
              // which restrains to use the same iterator loop logic
            }
        }

      if constexpr (contact_type == ContactType::particle_line ||
                    contact_type == ContactType::particle_point)
        {
          // Get current particle id
          unsigned int particle_id = pairs_in_contact_iterator->first;

          // Find if particle one is still a particle on this process
          auto particle_one_container = particle_container.find(particle_id);

          // Particle one is not on this process anymore, remove the container
          // and move on
          if (particle_one_container == particle_container.end())
            {
              pairs_in_contact_iterator =
                pairs_in_contact.erase(pairs_in_contact_iterator);
              continue;
            }

          // Get current particle id and particle iterator is updated
          adjacent_pairs_content->particle = particle_container.at(particle_id);
        }

      ++pairs_in_contact_iterator;
    }
}

// Particle container
template void
update_particle_container(
  DEM::dem_data_structures<2>::particle_index_iterator_map &particle_container,
  const Particles::ParticleHandler<2>                      *particle_handler);

template void
update_particle_container(
  DEM::dem_data_structures<3>::particle_index_iterator_map &particle_container,
  const Particles::ParticleHandler<3>                      *particle_handler);

// Local particle-particle contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::local_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Ghost particle-particle contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Local-local particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::local_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::local_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Local-ghost particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Ghost-local particle-particle periodic contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::adjacent_particle_pairs,
  ContactType::ghost_local_periodic_particle_particle>(
  DEM::dem_data_structures<2>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::adjacent_particle_pairs,
  ContactType::ghost_local_periodic_particle_particle>(
  DEM::dem_data_structures<3>::adjacent_particle_pairs &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-wall contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_wall>(
  DEM::dem_data_structures<2>::particle_wall_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_wall>(
  DEM::dem_data_structures<3>::particle_wall_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-floating wall contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  DEM::dem_data_structures<2>::particle_wall_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::particle_wall_in_contact,
  ContactType::particle_floating_wall>(
  DEM::dem_data_structures<3>::particle_wall_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-floating mesh contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::particle_triangle_cell_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  DEM::dem_data_structures<2>::particle_triangle_cell_from_mesh_in_contact
    &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::particle_triangle_cell_from_mesh_in_contact,
  ContactType::particle_floating_mesh>(
  DEM::dem_data_structures<3>::particle_triangle_cell_from_mesh_in_contact
    &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-point contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::particle_point_in_contact,
  ContactType::particle_point>(
  DEM::dem_data_structures<2>::particle_point_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::particle_point_in_contact,
  ContactType::particle_point>(
  DEM::dem_data_structures<3>::particle_point_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);

// Particle-line contact container
template void
update_contact_container_iterators<
  2,
  DEM::dem_data_structures<2>::particle_line_in_contact,
  ContactType::particle_line>(
  DEM::dem_data_structures<2>::particle_line_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<2>::particle_index_iterator_map
    &particle_container);

template void
update_contact_container_iterators<
  3,
  DEM::dem_data_structures<3>::particle_line_in_contact,
  ContactType::particle_line>(
  DEM::dem_data_structures<3>::particle_line_in_contact &pairs_in_contact,
  const DEM::dem_data_structures<3>::particle_index_iterator_map
    &particle_container);
