/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
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
 */

#include <dem/dem_contact_manager.h>

template <int dim>
void
DEMContactManager<dim>::execute_cell_neighbors_search(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
             periodic_boundaries_cells_information,
  const bool has_periodic_boundaries,
  const bool has_floating_mesh)
{
  // Find cell neighbors
  cell_neighbors_object.find_cell_neighbors(triangulation,
                                            cells_local_neighbor_list,
                                            cells_ghost_neighbor_list);

  // Find cell periodic neighbors
  if (has_periodic_boundaries)
    {
      cell_neighbors_object.find_cell_periodic_neighbors(
        triangulation,
        periodic_boundaries_cells_information,
        cells_local_periodic_neighbor_list,
        cells_ghost_periodic_neighbor_list,
        cells_ghost_local_periodic_neighbor_list);
    }

  // Get total (with repetition) neighbors list for floating mesh.
  if (has_floating_mesh)
    {
      cell_neighbors_object.find_full_cell_neighbors(triangulation,
                                                     total_neighbor_list);
    }
}

template <int dim>
void
DEMContactManager<dim>::update_contacts(const bool has_periodic_boundaries)
{
  // Update particle-particle contacts in local_adjacent_particles of fine
  // search step with local_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    typename DEM::dem_data_structures<dim>::particle_particle_candidates,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          local_contact_pair_candidates);

  // Update particle-particle contacts in global_adjacent_particles of fine
  // search step with global_contact_pair_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    typename DEM::dem_data_structures<dim>::particle_particle_candidates,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          ghost_contact_pair_candidates);

  if (has_periodic_boundaries)
    {
      // Update periodic particle-particle contacts in
      // local_periodic_adjacent_particles of fine search step with
      // local_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        typename DEM::dem_data_structures<dim>::particle_particle_candidates,
        ContactType::local_periodic_particle_particle>(
        local_periodic_adjacent_particles,
        local_contact_pair_periodic_candidates);

      // Update periodic particle-particle contacts in
      // ghost_periodic_adjacent_particles of fine search step with
      // ghost_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        typename DEM::dem_data_structures<dim>::particle_particle_candidates,
        ContactType::ghost_periodic_particle_particle>(
        ghost_periodic_adjacent_particles,
        ghost_contact_pair_periodic_candidates);

      // Update periodic particle-particle contacts in
      // ghost_local_periodic_adjacent_particles of fine search step with
      // ghost_local_contact_pair_periodic_candidates
      update_fine_search_candidates<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        typename DEM::dem_data_structures<dim>::particle_particle_candidates,
        ContactType::ghost_local_periodic_particle_particle>(
        ghost_local_periodic_adjacent_particles,
        ghost_local_contact_pair_periodic_candidates);
    }

  // Update particle-wall contacts in particle_wall_pairs_in_contact of fine
  // search step with particle_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates,
    ContactType::particle_wall>(particle_wall_in_contact,
                                particle_wall_candidates);

  // Update particle-floating wall contacts in particle_floating_wall_in_contact
  // of fine search step with particle_floating_wall_contact_candidates
  update_fine_search_candidates<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates,
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
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_candidates,
        ContactType::particle_floating_mesh>(
        particle_floating_mesh_in_contact[solid_counter],
        particle_floating_mesh_candidates[solid_counter]);
    }
}

template <int dim>
void
DEMContactManager<dim>::update_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler,
  const bool                             clear_contact_structures,
  const bool                             has_periodic_boundaries)
{
  // Update the iterators to local particles in a map of particles
  update_particle_container<dim>(particle_container, &particle_handler);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          particle_container,
                                          clear_contact_structures);

  // Update contact containers for ghost particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          particle_container,
                                          clear_contact_structures);

  if (has_periodic_boundaries)
    {
      // Update contact containers for local-local periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::local_periodic_particle_particle>(
        local_periodic_adjacent_particles,
        particle_container,
        clear_contact_structures);

      // Update contact containers for local-ghost periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::ghost_periodic_particle_particle>(
        ghost_periodic_adjacent_particles,
        particle_container,
        clear_contact_structures);

      // Update contact containers for ghost-local periodic particle-particle
      // pairs in contact
      update_contact_container_iterators<
        dim,
        typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
        ContactType::ghost_local_periodic_particle_particle>(
        ghost_local_periodic_adjacent_particles,
        particle_container,
        clear_contact_structures);
    }

  // Update contact containers for particle-wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
    ContactType::particle_wall>(particle_wall_in_contact, particle_container);

  // Update contact containers for particle-floating wall pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact,
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
        typename DEM::dem_data_structures<
          dim>::particle_floating_wall_from_mesh_in_contact,
        ContactType::particle_floating_mesh>(
        particle_floating_mesh_in_contact[solid_counter], particle_container);
    }

  // Update contact containers for particle-line pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_point_line_contact_info,
    ContactType::particle_point>(particle_lines_in_contact, particle_container);

  // Update contact containers for particle-point pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::particle_point_line_contact_info,
    ContactType::particle_point>(particle_points_in_contact,
                                 particle_container);
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_particle_broad_search(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const bool                               has_periodic_boundaries)
{
  particle_particle_broad_search_object.find_particle_particle_contact_pairs(
    particle_handler, *this);

  if (has_periodic_boundaries)
    {
      particle_particle_broad_search_object
        .find_particle_particle_periodic_contact_pairs(particle_handler, *this);
    }
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_particle_broad_search(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const DisableContacts<dim> &             disable_contacts_object,
  const bool                               has_periodic_boundaries)
{
  particle_particle_broad_search_object.find_particle_particle_contact_pairs(
    particle_handler, *this, disable_contacts_object);

  if (has_periodic_boundaries)
    {
      particle_particle_broad_search_object
        .find_particle_particle_periodic_contact_pairs(particle_handler,
                                                       *this,
                                                       disable_contacts_object);
    }
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_wall_broad_search(
  const Particles::ParticleHandler<dim> &particle_handler,
  BoundaryCellsInformation<dim> &        boundary_cell_object,
  const typename dem_data_structures<dim>::floating_mesh_information
                                                    floating_mesh_info,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
  const double                                      simulation_time,
  const bool                                        has_floating_mesh)
{
  // Particle-wall contact candidates
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cell_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_candidates);

  // Particle-floating wall contact pairs
  if (floating_walls.floating_walls_number > 0)
    {
      particle_wall_broad_search_object
        .find_particle_floating_wall_contact_pairs(
          boundary_cell_object.get_boundary_cells_with_floating_walls(),
          particle_handler,
          floating_walls,
          simulation_time,
          particle_floating_wall_candidates);
    }

  // Particle-floating mesh broad search
  if (has_floating_mesh)
    {
      particle_wall_broad_search_object.particle_floating_mesh_contact_search(
        floating_mesh_info,
        particle_handler,
        particle_floating_mesh_candidates,
        total_neighbor_list);
    }

  particle_point_candidates =
    particle_point_line_broad_search_object.find_particle_point_contact_pairs(
      particle_handler, boundary_cell_object.get_boundary_cells_with_points());

  if constexpr (dim == 3)
    {
      particle_line_candidates =
        particle_point_line_broad_search_object
          .find_particle_line_contact_pairs(
            particle_handler,
            boundary_cell_object.get_boundary_cells_with_lines());
    }
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_wall_broad_search(
  const Particles::ParticleHandler<dim> &particle_handler,
  BoundaryCellsInformation<dim> &        boundary_cell_object,
  const typename dem_data_structures<dim>::floating_mesh_information
                                                    floating_mesh_info,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
  const double                                      simulation_time,
  const DisableContacts<dim> &                      disable_contacts_object,
  const bool                                        has_floating_mesh)
{
  // Particle-wall contact candidates
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cell_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_candidates,
    disable_contacts_object);

  // Particle-floating wall contact pairs
  if (floating_walls.floating_walls_number > 0)
    {
      particle_wall_broad_search_object
        .find_particle_floating_wall_contact_pairs(
          boundary_cell_object.get_boundary_cells_with_floating_walls(),
          particle_handler,
          floating_walls,
          simulation_time,
          particle_floating_wall_candidates,
          disable_contacts_object);
    }

  // Particle-floating mesh broad search
  if (has_floating_mesh)
    {
      particle_wall_broad_search_object.particle_floating_mesh_contact_search(
        floating_mesh_info,
        particle_handler,
        particle_floating_mesh_candidates,
        total_neighbor_list,
        disable_contacts_object);
    }

  particle_point_candidates =
    particle_point_line_broad_search_object.find_particle_point_contact_pairs(
      particle_handler,
      boundary_cell_object.get_boundary_cells_with_points(),
      disable_contacts_object);

  if constexpr (dim == 3)
    {
      particle_line_candidates =
        particle_point_line_broad_search_object
          .find_particle_line_contact_pairs(
            particle_handler,
            boundary_cell_object.get_boundary_cells_with_lines(),
            disable_contacts_object);
    }
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_particle_fine_search(
  const double         neighborhood_threshold,
  const bool           has_periodic_boundaries,
  const Tensor<1, dim> periodic_offset)
{
  // Fine search for local particle-particle
  particle_particle_fine_search_object.particle_particle_fine_search(
    particle_container,
    local_adjacent_particles,
    local_contact_pair_candidates,
    neighborhood_threshold);

  // Fine search for ghost particle-particle
  particle_particle_fine_search_object.particle_particle_fine_search(
    particle_container,
    ghost_adjacent_particles,
    ghost_contact_pair_candidates,
    neighborhood_threshold);

  if (has_periodic_boundaries)
    {
      // Fine search for local-local periodic particle-particle
      particle_particle_fine_search_object.particle_particle_fine_search(
        particle_container,
        local_periodic_adjacent_particles,
        local_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);

      // Fine search for local-ghost periodic particle-particle
      particle_particle_fine_search_object.particle_particle_fine_search(
        particle_container,
        ghost_periodic_adjacent_particles,
        ghost_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);

      // Fine search for ghost-local periodic particle-particle
      particle_particle_fine_search_object.particle_particle_fine_search(
        particle_container,
        ghost_local_periodic_adjacent_particles,
        ghost_local_contact_pair_periodic_candidates,
        neighborhood_threshold,
        periodic_offset);
    }
}

template <int dim>
void
DEMContactManager<dim>::execute_particle_wall_fine_search(
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
  const double                                      simulation_time,
  const double                                      neighborhood_threshold,
  const bool                                        has_floating_mesh)
{
  // Particle - wall fine search
  particle_wall_fine_search_object.particle_wall_fine_search(
    particle_wall_candidates, particle_wall_in_contact);

  // Particle - floating wall fine search
  if (floating_walls.floating_walls_number > 0)
    {
      particle_wall_fine_search_object.particle_floating_wall_fine_search(
        particle_floating_wall_candidates,
        floating_walls,
        simulation_time,
        particle_floating_wall_in_contact);
    }

  // Particle - floating mesh fine search
  if (has_floating_mesh)
    {
      particle_wall_fine_search_object.particle_floating_mesh_fine_search(
        particle_floating_mesh_candidates, particle_floating_mesh_in_contact);
    }

  particle_points_in_contact =
    particle_point_line_fine_search_object.particle_point_fine_search(
      particle_point_candidates, neighborhood_threshold);

  if constexpr (dim == 3)
    {
      particle_lines_in_contact =
        particle_point_line_fine_search_object.particle_line_fine_search(
          particle_line_candidates, neighborhood_threshold);
    }
}

template class DEMContactManager<2>;
template class DEMContactManager<3>;
