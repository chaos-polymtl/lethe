// SPDX-FileCopyrightText: Copyright (c) 2022-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/dem_contact_manager.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/particle_wall_fine_search.h>
#include <dem/update_fine_search_candidates.h>
#include <dem/update_local_particle_containers.h>


using namespace DEM;

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::execute_cell_neighbors_search(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const typename dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information)
{
  // Get action manager
  auto *action_manager = DEMActionManager::get_action_manager();

  // Find cell neighbors
  find_cell_neighbors<dim>(triangulation,
                           cells_local_neighbor_list,
                           cells_ghost_neighbor_list);

  // Find cell periodic neighbors
  if (action_manager->check_periodic_boundaries_enabled())
    {
      find_cell_periodic_neighbors<dim>(
        triangulation,
        periodic_boundaries_cells_information,
        cells_local_periodic_neighbor_list,
        cells_ghost_periodic_neighbor_list,
        cells_ghost_local_periodic_neighbor_list);
    }

  // Get total (with repetition) neighbors list for floating mesh.
  if (action_manager->check_solid_objects_enabled())
    {
      find_full_cell_neighbors<dim>(triangulation, total_neighbor_list);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::update_contacts()
{
  // Update particle-particle contacts in local_adjacent_particles of fine
  // search step with local_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_structures<dim>::adjacent_particle_pairs,
    typename dem_data_structures<dim>::particle_particle_candidates,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          local_contact_pair_candidates);

  // Update particle-particle contacts in global_adjacent_particles of fine
  // search step with global_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_structures<dim>::adjacent_particle_pairs,
    typename dem_data_structures<dim>::particle_particle_candidates,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          ghost_contact_pair_candidates);

  if (DEMActionManager::get_action_manager()
        ->check_periodic_boundaries_enabled())
    {
      // Update periodic particle-particle contacts in
      // local_local_periodic_adjacent_particles of fine search step with
      // local_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        typename dem_data_structures<dim>::particle_particle_candidates,
        ContactType::local_periodic_particle_particle>(
        local_local_periodic_adjacent_particles,
        local_contact_pair_periodic_candidates);

      // Update periodic particle-particle contacts in
      // local_ghost_periodic_adjacent_particles of fine search step with
      // ghost_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        typename dem_data_structures<dim>::particle_particle_candidates,
        ContactType::ghost_periodic_particle_particle>(
        local_ghost_periodic_adjacent_particles,
        ghost_contact_pair_periodic_candidates);

      // Update periodic particle-particle contacts in
      // ghost_local_periodic_adjacent_particles of fine search step with
      // ghost_local_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        typename dem_data_structures<dim>::particle_particle_candidates,
        ContactType::ghost_local_periodic_particle_particle>(
        ghost_local_periodic_adjacent_particles,
        ghost_local_contact_pair_periodic_candidates);
    }

  // Update particle-wall contacts in particle_wall_pairs_in_contact of fine
  // search step with particle_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_structures<dim>::particle_wall_in_contact,
    typename dem_data_structures<dim>::particle_wall_candidates,
    ContactType::particle_wall>(particle_wall_in_contact,
                                particle_wall_candidates);

  // Update particle-floating wall contacts in particle_floating_wall_in_contact
  // of fine search step with particle_floating_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename dem_data_structures<dim>::particle_wall_in_contact,
    typename dem_data_structures<dim>::particle_floating_wall_candidates,
    ContactType::particle_floating_wall>(particle_floating_wall_in_contact,
                                         particle_floating_wall_candidates);

  // Update particle-floating mesh contacts in particle_floating_mesh_in_contact
  // of fine search step with particle_floating_mesh_contact_candidates
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_in_contact.size();
       ++solid_counter)
    {
      update_fine_search_candidates<
        dim,
        typename dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        typename dem_data_structures<
          dim>::particle_floating_wall_from_mesh_candidates,
        ContactType::particle_floating_mesh>(
        particle_floating_mesh_in_contact[solid_counter],
        particle_floating_mesh_candidates[solid_counter]);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::update_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  // Update the iterators to local particles in a map of particles
  update_particle_container<dim>(particle_container, &particle_handler);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          particle_container);

  // Update contact containers for ghost particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          particle_container);

  if (DEMActionManager::get_action_manager()
        ->check_periodic_boundaries_enabled())
    {
      // Update contact containers for local-local periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::local_periodic_particle_particle>(
        local_local_periodic_adjacent_particles, particle_container);

      // Update contact containers for local-ghost periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::ghost_periodic_particle_particle>(
        local_ghost_periodic_adjacent_particles, particle_container);

      // Update contact containers for ghost-local periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::ghost_local_periodic_particle_particle>(
        ghost_local_periodic_adjacent_particles, particle_container);
    }

  // Update contact containers for particle-wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::particle_wall_in_contact,
    ContactType::particle_wall>(particle_wall_in_contact, particle_container);

  // Update contact containers for particle-floating wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::particle_wall_in_contact,
    ContactType::particle_floating_wall>(particle_floating_wall_in_contact,
                                         particle_container);

  // Update contact containers for particle-floating mesh pairs in contact for
  // every solid objects of mesh
  for (unsigned int solid_counter = 0;
       solid_counter < particle_floating_mesh_in_contact.size();
       ++solid_counter)
    {
      update_contact_container_iterators<
        dim,
        typename dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        ContactType::particle_floating_mesh>(
        particle_floating_mesh_in_contact[solid_counter], particle_container);
    }

  // Update contact containers for particle-line pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::particle_line_in_contact,
    ContactType::particle_line>(particle_lines_in_contact, particle_container);

  // Update contact containers for particle-point pairs in contact
  update_contact_container_iterators<
    dim,
    typename dem_data_structures<dim>::particle_point_in_contact,
    ContactType::particle_point>(particle_points_in_contact,
                                 particle_container);
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::execute_particle_particle_broad_search(
  dealii::Particles::ParticleHandler<dim>            &particle_handler,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  auto *action_manager = DEMActionManager::get_action_manager();

  // Check if sparse contacts are enabled to use proper broad search functions
  // The first broad search is the default one for sparse contacts
  if (action_manager->use_default_broad_search_functions())
    {
      find_particle_particle_contact_pairs<dim>(particle_handler,
                                                cells_local_neighbor_list,
                                                cells_ghost_neighbor_list,
                                                local_contact_pair_candidates,
                                                ghost_contact_pair_candidates);

      if (action_manager->check_periodic_boundaries_enabled())
        {
          find_particle_particle_periodic_contact_pairs<dim>(
            particle_handler,
            cells_local_periodic_neighbor_list,
            cells_ghost_periodic_neighbor_list,
            cells_ghost_local_periodic_neighbor_list,
            local_contact_pair_periodic_candidates,
            ghost_contact_pair_periodic_candidates,
            ghost_local_contact_pair_periodic_candidates);
        }
    }
  else
    {
      find_particle_particle_contact_pairs<dim>(particle_handler,
                                                cells_local_neighbor_list,
                                                cells_ghost_neighbor_list,
                                                local_contact_pair_candidates,
                                                ghost_contact_pair_candidates,
                                                sparse_contacts_object);

      if (action_manager->check_periodic_boundaries_enabled())
        {
          find_particle_particle_periodic_contact_pairs<dim>(
            particle_handler,
            cells_local_periodic_neighbor_list,
            cells_ghost_periodic_neighbor_list,
            cells_ghost_local_periodic_neighbor_list,
            local_contact_pair_periodic_candidates,
            ghost_contact_pair_periodic_candidates,
            ghost_local_contact_pair_periodic_candidates,
            sparse_contacts_object);
        }
    }
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::execute_particle_wall_broad_search(
  const Particles::ParticleHandler<dim> &particle_handler,
  BoundaryCellsInformation<dim>         &boundary_cell_object,
  const typename dem_data_structures<dim>::solid_surfaces_mesh_information
                                                      solid_surfaces_mesh_info,
  const Parameters::Lagrangian::FloatingWalls<dim>   &floating_walls,
  const double                                        simulation_time,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object)
{
  auto *action_manager = DEMActionManager::get_action_manager();

  if (action_manager->use_default_broad_search_functions())
    {
      // Particle-wall contact candidates
      find_particle_wall_contact_pairs<dim>(
        boundary_cell_object.get_boundary_cells_information(),
        particle_handler,
        particle_wall_candidates);

      // Particle-floating wall contact pairs
      if (floating_walls.floating_walls_number > 0)
        {
          find_particle_floating_wall_contact_pairs<dim>(
            boundary_cell_object.get_boundary_cells_with_floating_walls(),
            particle_handler,
            floating_walls,
            simulation_time,
            particle_floating_wall_candidates);
        }

      // Particle-floating mesh broad search
      if (action_manager->check_solid_objects_enabled())
        {
          particle_solid_surfaces_contact_search<dim>(
            solid_surfaces_mesh_info,
            particle_handler,
            particle_floating_mesh_candidates,
            total_neighbor_list);
        }

      find_particle_point_contact_pairs<dim>(
        particle_handler,
        boundary_cell_object.get_boundary_cells_with_points(),
        particle_point_candidates);

      if constexpr (dim == 3)
        {
          find_particle_line_contact_pairs<dim>(
            particle_handler,
            boundary_cell_object.get_boundary_cells_with_lines(),
            particle_line_candidates);
        }
    }
  else
    {
      // Particle-wall contact candidates
      find_particle_wall_contact_pairs<dim>(
        boundary_cell_object.get_boundary_cells_information(),
        particle_handler,
        particle_wall_candidates,
        sparse_contacts_object);

      // Particle-floating wall contact pairs
      if (floating_walls.floating_walls_number > 0)
        {
          find_particle_floating_wall_contact_pairs<dim>(
            boundary_cell_object.get_boundary_cells_with_floating_walls(),
            particle_handler,
            floating_walls,
            simulation_time,
            particle_floating_wall_candidates,
            sparse_contacts_object);
        }

      // Particle-floating mesh broad search
      if (action_manager->check_solid_objects_enabled())
        {
          particle_solid_surfaces_contact_search<dim>(
            solid_surfaces_mesh_info,
            particle_handler,
            particle_floating_mesh_candidates,
            total_neighbor_list,
            sparse_contacts_object);
        }

      find_particle_point_contact_pairs<dim>(
        particle_handler,
        boundary_cell_object.get_boundary_cells_with_points(),
        particle_point_candidates,
        sparse_contacts_object);

      if constexpr (dim == 3)
        {
          find_particle_line_contact_pairs<dim>(
            particle_handler,
            boundary_cell_object.get_boundary_cells_with_lines(),
            particle_line_candidates,
            sparse_contacts_object);
        }
    }
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::execute_particle_particle_fine_search(
  const double neighborhood_threshold)
{
  // Fine search for local particle-particle
  particle_particle_fine_search<dim>(particle_container,
                                     local_adjacent_particles,
                                     local_contact_pair_candidates,
                                     neighborhood_threshold);

  // Fine search for ghost particle-particle
  particle_particle_fine_search<dim>(particle_container,
                                     ghost_adjacent_particles,
                                     ghost_contact_pair_candidates,
                                     neighborhood_threshold);

  if (DEMActionManager::get_action_manager()
        ->check_periodic_boundaries_enabled())
    {
      // Fine search for local-local periodic particle-particle
      particle_particle_fine_search<dim>(
        particle_container,
        local_local_periodic_adjacent_particles,
        local_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);

      // Fine search for local-ghost periodic particle-particle
      particle_particle_fine_search<dim>(
        particle_container,
        local_ghost_periodic_adjacent_particles,
        ghost_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);

      // Fine search for ghost-local periodic particle-particle
      particle_particle_fine_search<dim>(
        particle_container,
        ghost_local_periodic_adjacent_particles,
        ghost_local_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMContactManager<dim, PropertiesIndex>::execute_particle_wall_fine_search(
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
  const double                                      simulation_time,
  const double                                      neighborhood_threshold)
{
  // Particle - wall fine search
  particle_wall_fine_search<dim>(particle_wall_candidates,
                                 particle_wall_in_contact);

  // Particle - floating wall fine search
  if (floating_walls.floating_walls_number > 0)
    {
      particle_floating_wall_fine_search<dim>(
        particle_floating_wall_candidates,
        floating_walls,
        simulation_time,
        particle_floating_wall_in_contact);
    }

  // Particle - floating mesh fine search
  if (DEMActionManager::get_action_manager()->check_solid_objects_enabled())
    {
      particle_floating_mesh_fine_search<dim>(
        particle_floating_mesh_candidates, particle_floating_mesh_in_contact);
    }

  particle_point_fine_search<dim, PropertiesIndex>(particle_point_candidates,
                                                   neighborhood_threshold,
                                                   particle_points_in_contact);

  if constexpr (dim == 3)
    {
      particle_line_fine_search<dim, PropertiesIndex>(
        particle_line_candidates,
        neighborhood_threshold,
        particle_lines_in_contact);
    }
}

template class DEMContactManager<2, DEM::DEMProperties::PropertiesIndex>;
template class DEMContactManager<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class DEMContactManager<2, DEM::DEMMPProperties::PropertiesIndex>;
template class DEMContactManager<3, DEM::DEMProperties::PropertiesIndex>;
template class DEMContactManager<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class DEMContactManager<3, DEM::DEMMPProperties::PropertiesIndex>;
