// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/find_cell_neighbors.h>
#include <dem/particle_ray_tracing.h>
#include <dem/read_mesh.h>

#include <sys/stat.h>
#include <sstream>


template <int dim, typename PropertiesIndex>
ParticleRayTracing<dim, PropertiesIndex>::ParticleRayTracing(
  DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , parameters(dem_parameters)
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , background_dh(triangulation)
  , photon_handler(triangulation, mapping, 3)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
//, photon_displacement_tensor(GridTools::minimal_cell_diameter(triangulation))
{
  Tensor<1, dim> temp;
  photon_displacement_tensor =
    temp * GridTools::minimal_cell_diameter(triangulation);
}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::setup_parameters()
{
  // Print simulation starting information
  pcout << std::endl;
  std::stringstream ss;
  ss << "Running on " << n_mpi_processes << " rank(s)";
  announce_string(pcout, ss.str(), '*');

  // Check if the output directory exists
  std::string output_dir_name = parameters.simulation_control.output_folder;
  struct stat buffer;

  // If output directory does not exist, create it
  if (this_mpi_process == 0)
    {
      if (stat(output_dir_name.c_str(), &buffer) != 0)
        {
          create_output_folder(output_dir_name);
        }
    }

}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::find_locally_own_cells_near_particles(
  typename dem_data_structures<dim>::cell_set local_cells_near_particles)
{
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;

  // Build the neighboring cell lists
  find_cell_neighbors<dim, true>(triangulation,
                                 cells_local_neighbor_list,
                                 cells_ghost_neighbor_list);

  // Loop over every vector in the local-cell -> local-neighbor-cell container
  for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_local_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The first iterator is the main_cell itself.
      auto main_cell_iterator = cell_neighbor_list_iterator->begin();

      // Loop over the local neighboring cells.
      for (auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
           cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_list_iterator)
        {
          // If the neighboring cell contains at least one particle we add the
          // main cell in the set
          if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) !=
              0)
            {
              local_cells_near_particles.insert(*main_cell_iterator);
              continue;
            }
        }
    }

  // Loop over every vector in the local-cell -> ghost-neighbor-cell container
  // With this loop, we flag local cell with ghost neighbors containing
  // particles.
  for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The first iterator is the main_cell itself.
      auto main_cell_iterator = cell_neighbor_list_iterator->begin();

      // Check if the cell is already in the set from the previous loop.
      auto candidates_container_it =
        local_cells_near_particles.find(*main_cell_iterator);

      if (candidates_container_it == local_cells_near_particles.end())
        continue;

      // Loop over the local neighboring cells.
      for (auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
           cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_list_iterator)
        {
          if (particle_handler.n_particles_in_cell(*cell_neighbor_iterator) !=
              0)
            {
              local_cells_near_particles.insert(*main_cell_iterator);
              continue;
            }
        }
    }
}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::insert_photon()
{
  // Vector taken from the dem_parameters (unit vectors)
  const Tensor<1, dim> dir_1;
  const Tensor<1, dim> dir_2;

  // Point taken from the dem_parameters
  const Point<dim> origin;

  // Delta in each direction between each photon
  const double delta_1 = 1.0;
  const double delta_2 = 1.1;

  // Double taken from the dem_parameters (number of photon per direction)
  const unsigned int n_1 = 10;
  const unsigned int n_2 = 11;

  unsigned int n_total_particles_to_insert = n_1 * n_2;
  unsigned int n_particle_to_insert_on_this_proc, start_id;

  const unsigned int base      = n_total_particles_to_insert / n_mpi_processes;
  const unsigned int remainder = n_total_particles_to_insert % n_mpi_processes;

  if (this_mpi_process < remainder)
    {
      n_particle_to_insert_on_this_proc = base + 1;
      start_id                          = this_mpi_process * (base + 1);
    }
  else
    {
      n_particle_to_insert_on_this_proc = base;
      start_id = remainder * (base + 1) + (this_mpi_process - remainder) * base;
    }

  // Insertion point
  std::vector<Point<dim>> insertion_points_on_proc;

  // Photon properties are the insertion point position. This is used to compute
  // the closest insertion point with the origin of the photon/ray.
  std::vector<std::vector<double>> photon_properties;

  insertion_points_on_proc.reserve(n_particle_to_insert_on_this_proc);
  photon_properties.reserve(n_particle_to_insert_on_this_proc);
  Point<dim>          insertion_point_it;
  std::vector<double> properties_of_one_photon;

  // For now, insertion is only carried by processor 0.
  for (unsigned int k = 0; k < n_particle_to_insert_on_this_proc; ++k)
    {
      unsigned int global_index = start_id + k;
      unsigned int i            = global_index / n_2;
      unsigned int j            = global_index % n_2;

      insertion_point_it = origin + i * delta_1 * dir_1 + j * delta_2 * dir_2;

      insertion_points_on_proc.push_back(insertion_point_it);

      properties_of_one_photon = [&]() {
        if constexpr (dim == 2)
          {
            return std::vector<double>{insertion_point_it(0), insertion_point_it(1)};
          }
        else if constexpr (dim == 3)
          {
            return std::vector<double>{insertion_point_it(0),
                    insertion_point_it(1),
                    insertion_point_it(2)};
          }
      }();

      photon_properties.push_back(properties_of_one_photon);
    }

  MPI_Comm           communicator = triangulation.get_communicator();
  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);
  // Obtaining global bounding boxes
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(communicator, my_bounding_box);

  pcout << "Inserting photons" << std::endl;

  // Insert the photons using the insertion points and properties
  photon_handler.insert_global_particles(insertion_points_on_proc,
                                         global_bounding_boxes,
                                         photon_properties);
}

template <int dim, typename PropertiesIndex>
void
ParticleRayTracing<dim, PropertiesIndex>::solve()
{
  // Triangulation
  // Insert particle
  // Insert photons

  // While number of photon bigger than 0

  //  // Pseudo integrate

  //  // Store contact intersection point

  // Output

}

template class ParticleRayTracing<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleRayTracing<2, DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleRayTracing<3, DEM::DEMMPProperties::PropertiesIndex>;
