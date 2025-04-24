// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/dem_action_manager.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_ray_tracing.h>


template <int dim, typename PropertiesIndex>
ParticleRayTracing<dim, PropertiesIndex>::ParticleRayTracing(
  DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , background_dh(triangulation)
  , photon_handler(triangulation, mapping, 1)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
{}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::find_locally_own_cells_with_particles(
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_with_particle_and_neighbours)
{
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;

  find_cell_neighbors<dim, true>(triangulation,
                                 cells_local_neighbor_list,
                                 cells_ghost_neighbor_list);

  // Loop over every vector in the local-cell -> local-neighbour-cell container
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      auto main_cell_iterator = cell_neighbor_list_iterator->begin();

      // Loop over the local neighboring cells.
      // The first iterator is the main_cell itself.
      for (auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
           cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_list_iterator)
        {
          if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) !=
              0)
            {
              local_and_ghost_cells_with_particles_and_neighbors.insert(
                *main_cell_iterator);
              continue;
            }
        }
    }

  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      auto main_cell_iterator = cell_neighbor_list_iterator->begin();

      auto candidates_container_it =
        local_and_ghost_cells_with_particles_and_neighbors.find(
          *main_cell_iterator);

      // Check if the cell is already in the set from the previous loop.
      if (candidates_container_it ==
          local_and_ghost_cells_with_particles_and_neighbors.end())
        continue;

      // Loop over the local neighboring cells.
      // The first iterator is the main_cell itself.
      for (auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
           cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_list_iterator)
        {
          if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) !=
              0)
            {
              local_and_ghost_cells_with_particles_and_neighbors.insert(
                *main_cell_iterator);
              continue;
            }
        }
    }
}



template class ParticleRayTracing<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<2, DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMMPProperties::PropertiesIndex>;
