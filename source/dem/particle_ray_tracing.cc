// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_ray_tracing.h>


template <int dim, typename PropertiesIndex>
ParticleRayTracing<dim, PropertiesIndex>::ParticleRayTracing(
  ParticleRayTracingParameters<dim> parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , photon_handler(triangulation, mapping, 1)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
{}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::find_locally_own_cells_with_particles(
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_with_particle_and_neighbours)
{
  /* When we will compute the intersection with
  // Get action manager
  auto *action_manager = DEMActionManager::get_action_manager();

 // Get total (with repetition) neighbors list for floating mesh.
 if (action_manager->check_solid_objects_enabled())
   {
     find_full_cell_neighbors<dim>(triangulation, total_neighbor_list);
   }
 */
  find_cell_neighbors<dim>(triangulation,
                           cells_local_neighbor_list,
                           cells_ghost_neighbor_list);

  // Loop over every vector in the local-cell -> local-neighbour-cell container
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
          // We add the main cell to the
          local_and_ghost_cells_with_particles_and_neighbors.insert(*cell_neighbor_iterator);

          // Going through neighbor cells of the main cell
          ++cell_neighbor_iterator;
          for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
               ++cell_neighbor_iterator)
            {
              // We check if the neighbors cells are already in the set.
              // if not, we add it.
              auto search_iterator =
                local_and_ghost_cells_with_particles_and_neighbors.find(
                  *cell_neighbor_iterator);

              if (search_iterator == local_and_ghost_cells_with_particles_and_neighbors.end())
                local_and_ghost_cells_with_particles_and_neighbors.insert(*cell_neighbor_iterator);
            }
        }
    }
}



template class ParticleRayTracing<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<2, DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMMPProperties::PropertiesIndex>;
