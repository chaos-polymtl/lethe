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

#include <dem/dem_container_manager.h>

template <int dim>
void
DEMContainerManager<dim>::execute_particle_particle_board_search(
  dealii::Particles::ParticleHandler<dim> &particle_handler)
{
  // Clearing local & ghost contact pair candidates
  clear_contact_pair_candidates();

  // First we handle the local-local candidate pairs

  // Looping over the potential cells which may contain particles.
  // This includes the cell itself as well as the neighbouring cells that
  // were identified.
  // cell_neighbor_list_iterator is [cell_it, neighbor_0_it, neighbor_1_it, ...]
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      // Check to see if the main cell has any particles
      if (particles_exist_in_main_cell)
        {
          // Find local-local collision pairs in the main cell, 1st particle
          // iterator is skipped since the main particle will not be
          // considered as collision partner with itself
          for (auto particle_in_main_cell = particles_in_main_cell.begin();
               particle_in_main_cell != particles_in_main_cell.end();
               ++particle_in_main_cell)
            {
              store_candidates(
                std::next(particle_in_main_cell, 1),
                particles_in_main_cell,
                local_contact_pair_candidates[particle_in_main_cell->get_id()]);
            }

          // Going through neighbor cells of the main cell
          ++cell_neighbor_iterator;
          for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
               ++cell_neighbor_iterator)
            {
              // Defining iterator on local particles in the neighbor cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_neighbor_cell =
                  particle_handler.particles_in_cell(*cell_neighbor_iterator);

              // Capturing particle pairs, the first particle in the main
              // cell and the second particle in the neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates(particles_in_neighbor_cell.begin(),
                                   particles_in_neighbor_cell,
                                   local_contact_pair_candidates
                                     [particle_in_main_cell->get_id()]);
                }
            }
        }
    }

  // Now we go through the local-ghost pairs (the first iterator shows a local
  // particles, and the second a ghost particle)

  // Looping over cells_ghost_neighbor_list
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();

      // Particles in the main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_main_cell =
          particle_handler.particles_in_cell(*cell_neighbor_iterator);

      const bool particles_exist_in_main_cell = !particles_in_main_cell.empty();

      if (particles_exist_in_main_cell)
        {
          // Going through ghost neighbor cells of the main cell
          ++cell_neighbor_iterator;

          for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
               ++cell_neighbor_iterator)
            {
              // Defining iterator on ghost particles in the neighbor cells
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particles_in_neighbor_cell =
                  particle_handler.particles_in_cell(*cell_neighbor_iterator);

              // Capturing particle pairs, the first particle (local) in
              // the main cell and the second particle (ghost) in the
              // neighbor cells
              for (auto particle_in_main_cell = particles_in_main_cell.begin();
                   particle_in_main_cell != particles_in_main_cell.end();
                   ++particle_in_main_cell)
                {
                  store_candidates(particles_in_neighbor_cell.begin(),
                                   particles_in_neighbor_cell,
                                   ghost_contact_pair_candidates
                                     [particle_in_main_cell->get_id()]);
                }
            }
        }
    }
}

template <int dim>
void
DEMContainerManager<dim>::localize_contacts()
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
DEMContainerManager<dim>::locate_local_particles_in_cells(
  const Particles::ParticleHandler<dim> &particle_handler)
{
  // Update the iterators to local particles in a map of particles
  update_particle_container<dim>(particle_container, &particle_handler);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::local_particle_particle>(local_adjacent_particles,
                                          particle_container);

  // Update contact containers for local particle-particle pairs in contact
  update_contact_container_iterators<
    dim,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs,
    ContactType::ghost_particle_particle>(ghost_adjacent_particles,
                                          particle_container);

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

template class DEMContainerManager<2>;
template class DEMContainerManager<3>;
