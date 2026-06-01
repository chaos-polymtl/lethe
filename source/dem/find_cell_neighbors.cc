// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include "core/lethe_grid_tools.h"

#include <dem/find_cell_neighbors.h>

#include <deal.II/grid/grid_tools.h>

using namespace DEM;

template <int dim, bool reciprocal>
void
find_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list)
{
  // The output vectors of the function are cells_local_neighbor_list and
  // cells_ghost_neighbor_list. The first contains all the local neighbors cells
  // of all local cells; while the second contains all the ghost cells of all
  // local cells. They are two vectors with the size of the number of active
  // cells. The first elements of all vectors are the main cells
  typename dem_data_structures<dim>::cell_vector local_neighbor_vector;
  typename dem_data_structures<dim>::cell_vector ghost_neighbor_vector;

  // This vector is used to avoid repetition of adjacent cells. For instance if
  // cell B is recognized as the neighbor of cell A, cell A will not be added to
  // the neighbor list of cell B again. This is done using the total_cell_list
  // vector
  typename dem_data_structures<dim>::cell_set total_cell_list;

  // For each cell, the cell vertices are found and used to find adjacent cells.
  // The reason is to find the cells located on the corners of the main cell.
  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  // Looping over cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
          // The first element of each vector is the cell itself.
          local_neighbor_vector.push_back(cell);

          if constexpr (!reciprocal)
            total_cell_list.insert(cell);

          for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
            {
              for (const auto &neighbor : v_to_c[cell->vertex_index(vertex)])
                {
                  if (neighbor->is_locally_owned())
                    {
                      auto search_iterator = total_cell_list.find(neighbor);
                      auto local_search_iterator =
                        std::find(local_neighbor_vector.begin(),
                                  local_neighbor_vector.end(),
                                  neighbor);

                      // If the cell (neighbor) is a local cell and not present
                      // in the total_cell_list vector, it will be added as the
                      // neighbor of the main cell
                      // ("cell") and also to the total_cell_list to avoid
                      // repetition for next cells.
                      if (search_iterator == total_cell_list.end() &&
                          local_search_iterator == local_neighbor_vector.end())
                        {
                          local_neighbor_vector.push_back(neighbor);
                        }

                      // If the neighbor cell is a ghost, it should be added in
                      // the ghost_neighbor_vector container
                    }
                  else if (neighbor->is_ghost())
                    {
                      auto ghost_search_iterator =
                        std::find(ghost_neighbor_vector.begin(),
                                  ghost_neighbor_vector.end(),
                                  neighbor);
                      if (ghost_search_iterator == ghost_neighbor_vector.end())
                        {
                          if (ghost_neighbor_vector.empty())
                            {
                              ghost_neighbor_vector.push_back(cell);
                            }

                          ghost_neighbor_vector.push_back(neighbor);
                        }
                    }
                }
            }
        }
      if (!local_neighbor_vector.empty())
        cells_local_neighbor_list.push_back(local_neighbor_vector);
      if (!ghost_neighbor_vector.empty())
        cells_ghost_neighbor_list.push_back(ghost_neighbor_vector);

      local_neighbor_vector.clear();
      ghost_neighbor_vector.clear();
    }
}

template <int dim>
void
find_cell_periodic_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  const typename DEM::dem_data_structures<dim>::cell_touch_boundary_id
                                           &cell_to_pbc_mesh_id_set,
  const PeriodicBoundariesManipulator<dim> &periodic_boundaries_object,
  std::vector<typename DEM::dem_data_structures<dim>::cells_neighbor_list>
    &cells_local_periodic_neighbor_lists,
  std::vector<typename DEM::dem_data_structures<dim>::cells_neighbor_list>
    &cells_ghost_periodic_neighbor_lists,
  std::vector<typename DEM::dem_data_structures<dim>::cells_neighbor_list>
    &cells_ghost_local_periodic_neighbor_lists)
{
  std::map<std::uint8_t, typename dem_data_structures<dim>::cell_vector>
    local_periodic_neighbor_vectors, ghost_periodic_neighbor_vectors,
    ghost_local_periodic_neighbor_vectors;

  typename dem_data_structures<dim>::cell_set total_cell_list;
  typename dem_data_structures<dim>::cell_set total_ghost_cell_list;

  // For each cell, the cell vertices are found and used to find adjacent cells.
  // The reason is to find the cells located on the corners of the main cell.
  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  // A map of coinciding vertices labeled by an arbitrary element from them
  std::map<unsigned int, std::vector<unsigned int>> coinciding_vertex_groups;

  // Map of a vertex to the label of a group of coinciding vertices
  std::map<unsigned int, unsigned int> vertex_to_coinciding_vertex_group;

  // Collect for a given triangulation all locally relevant vertices that
  // coincide due to periodicity.
  GridTools::collect_coinciding_vertices(triangulation,
                                         coinciding_vertex_groups,
                                         vertex_to_coinciding_vertex_group);

  // Keep track of cells already processed to avoid redundant corner handling,
  // since periodic_boundaries_cells_information is an unordered_multimap
  std::set<typename Triangulation<dim>::active_cell_iterator> processed_cells;

  // Looping on every periodic main cell
  for (const auto &[main_cell, main_cell_touching_boundaries] :
       cell_to_pbc_mesh_id_set)
    {
      // Skip if we already mapped this primary cell.
      // NOTE: Since we are looping of the cell_to_pbc_mesh_id_set container, I
      // think this if is no longer required since we are looping on
      // cell_to_pbc_mesh_id_set. Each cell has one entry inside this container,
      // which wasn't the case when using the periodic info struc, since a cell
      // touching two periodic boundaries has two structs.
      if (processed_cells.contains(main_cell))
        continue;

      processed_cells.insert(main_cell);

      // If the main cell is owned by the processor
      if (main_cell->is_locally_owned())
        {
          // Store every main cell
          total_cell_list.insert(main_cell);

          // Empty list of periodic cell neighbors
          typename dem_data_structures<dim>::cell_vector
            periodic_neighbor_vector;

          // Get the periodic neighbor(s) of the main cell
          LetheGridTools::get_periodic_neighbor_list<dim>(
            main_cell,
            coinciding_vertex_groups,
            vertex_to_coinciding_vertex_group,
            v_to_c,
            periodic_neighbor_vector);

          // Loop on every periodic neighboring cell, determine in which
          // container it should go
          for (const auto &periodic_neighbor : periodic_neighbor_vector)
            {
              // Find the boundaries that are touching the current neighboring
              // cell.
              const std::set<types::boundary_id>
                &neighboring_cell_touching_boundaries =
                  cell_to_pbc_mesh_id_set.at(periodic_neighbor);

              // Find which boundary (or boundaries) is shared between the main
              // cell and its neighboring cell
              std::set<types::boundary_id> combination_of_shared_boundaries =
                periodic_boundaries_object.find_shared_periodic_boundaries(
                  main_cell_touching_boundaries,
                  neighboring_cell_touching_boundaries);

              // From the shared periodic boundaries 0 between the main and
              // neighboring cell, we find in which container the neighboring
              // cell should be inserted in.
              std::uint8_t container_index =
                periodic_boundaries_object.get_container_index(
                  combination_of_shared_boundaries);

              // If the neighboring cell is locally owned
              if (periodic_neighbor->is_locally_owned())
                {
                  // Check if the neighbor cell has been processed as a main
                  // cell.
                  auto search_iterator =
                    total_cell_list.find(periodic_neighbor);

                  // If the neighboring cell is a local cell and has still not
                  // been treated as a main cell, it will be added as to the
                  // neighboring cells of the current main cell.
                  if (search_iterator == total_cell_list.end())
                    {
                      // Check if
                      // local_periodic_neighbor_vectors[container_index] exist,
                      // if not, create it and add the main cell inside, since
                      // the main cell is always the first iterator of the
                      // vector.
                      if (!local_periodic_neighbor_vectors.contains(
                            container_index))
                        local_periodic_neighbor_vectors.insert(
                          {container_index, {main_cell}});

                      // Add the neighboring cell to the vectors
                      local_periodic_neighbor_vectors[container_index]
                        .push_back(periodic_neighbor);
                    }
                }
              else if (periodic_neighbor->is_ghost())
                {
                  // If the cell neighbor is a ghost, it should be added in
                  // ghost_periodic_neighbor_vectors[container_index]

                  // Check if
                  // ghost_periodic_neighbor_vectors[container_index] exist,
                  // if not, create it and add the main cell inside, since
                  // the main cell is always the first iterator of the vector.
                  if (!ghost_periodic_neighbor_vectors.contains(
                        container_index))
                    ghost_periodic_neighbor_vectors.insert(
                      {container_index, {main_cell}});

                  // Add the neighboring cell to the vectors
                  ghost_periodic_neighbor_vectors[container_index].push_back(
                    periodic_neighbor);
                }
            }

          // Now that every neighboring cell of the main cell has been treated,
          // we store the cell_vectors in the appropriate containers
          if (!local_periodic_neighbor_vectors.empty())
            {
              for (const auto &[key, value] : local_periodic_neighbor_vectors)
                cells_local_periodic_neighbor_lists[key].push_back(value);
              local_periodic_neighbor_vectors.clear();
            }

          if (!ghost_periodic_neighbor_vectors.empty())
            {
              for (const auto &[key, value] : ghost_periodic_neighbor_vectors)
                cells_ghost_periodic_neighbor_lists[key].push_back(value);
              ghost_periodic_neighbor_vectors.clear();
            }
        }
      else if (main_cell->is_ghost())
        {
          // Since periodic cells are mapped on one side only (cells on pb 0
          // with cells on pb 1), we need a 3rd type of container for
          // ghost-local contacts for force calculation, where the local cell is
          // on the pb 1 side. Here we store local neighbors of ghost cells.

          // Empty list of periodic cell neighbor
          typename dem_data_structures<dim>::cell_vector periodic_neighbor_list;

          // Get the periodic neighbor of the cell
          LetheGridTools::get_periodic_neighbor_list<dim>(
            main_cell,
            coinciding_vertex_groups,
            vertex_to_coinciding_vertex_group,
            v_to_c,
            periodic_neighbor_list);

          for (const auto &periodic_neighbor : periodic_neighbor_list)
            {
              if (periodic_neighbor->is_locally_owned())
                {
                  // Find the boundaries that are touching the current
                  // neighboring cell.
                  const std::set<types::boundary_id>
                    neighboring_cell_touching_boundaries =
                      cell_to_pbc_mesh_id_set.at(periodic_neighbor);

                  // Find which boundary (or boundaries) is shared between
                  // the main cell and its neighboring cell
                  std::set<types::boundary_id>
                    combination_of_shared_boundaries =
                      periodic_boundaries_object
                        .find_shared_periodic_boundaries(
                          main_cell_touching_boundaries,
                          neighboring_cell_touching_boundaries);

                  // From the shared boundaries, we find in which container
                  // the neighboring cell should be inserted in.
                  std::uint8_t container_index =
                    periodic_boundaries_object.get_container_index(
                      combination_of_shared_boundaries);

                  if (!ghost_periodic_neighbor_vectors.contains(
                        container_index))
                    ghost_local_periodic_neighbor_vectors[container_index]
                      .push_back(main_cell);

                  ghost_local_periodic_neighbor_vectors[container_index]
                    .push_back(periodic_neighbor);
                }
            }

          // Loop on the ghost_local_periodic_neighbor_vector key and values.
          // Insert the value inside the
          // ells_ghost_local_periodic_neighbor_lists[key]
          if (!ghost_local_periodic_neighbor_vectors.empty())
            {
              for (const auto &[key, value] :
                   ghost_local_periodic_neighbor_vectors)
                {
                  cells_ghost_local_periodic_neighbor_lists[key].push_back(
                    value);
                }
            }
        }
    }
}

// This function finds the full neighbor list (with repetition) of all the
// active cells in the triangulation. Because particle-floating mesh
// contacts need all the particles located in ALL the neighbor cells of
// the main cell to search for possible collisions with the floating mesh.
template <int dim>
void
find_full_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename dem_data_structures<dim>::cells_total_neighbor_list
    &cells_total_neighbor_list)
{
  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  // Looping over cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
          std::vector<typename Triangulation<dim>::active_cell_iterator>
            full_neighbor_vector(0);

          full_neighbor_vector.push_back(cell);

          for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
            {
              for (const auto &neighbor : v_to_c[cell->vertex_index(vertex)])
                {
                  if (neighbor->is_locally_owned())
                    {
                      auto search_iterator =
                        std::find(full_neighbor_vector.begin(),
                                  full_neighbor_vector.end(),
                                  neighbor);

                      if (search_iterator == full_neighbor_vector.end())
                        {
                          full_neighbor_vector.push_back(neighbor);
                        }
                    }
                }
            }

          cells_total_neighbor_list.insert(
            {cell->global_active_cell_index(), full_neighbor_vector});
        }
    }
}

// template <int dim>
// void
// get_periodic_neighbor_list(
//   const typename Triangulation<dim>::active_cell_iterator &cell,
//   const std::map<unsigned int, std::vector<unsigned int>>
//                                              &coinciding_vertex_groups,
//   const std::map<unsigned int, unsigned int>
//   &vertex_to_coinciding_vertex_group, const std::vector<std::set<typename
//   Triangulation<dim>::active_cell_iterator>>
//                                                  &v_to_c,
//   typename dem_data_structures<dim>::cell_vector &periodic_neighbor_list)
// {
//   std::set<typename Triangulation<dim>::active_cell_iterator>
//   unique_neighbors;
//
//   // Loop over all vertices of the cell
//   for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
//     {
//       const unsigned int vertex_id = cell->vertex_index(vertex);
//
//       auto group_it = vertex_to_coinciding_vertex_group.find(vertex_id);
//
//       if (group_it == vertex_to_coinciding_vertex_group.end())
//         continue;
//
//       const unsigned int coinciding_vertex_key = group_it->second;
//
//       for (const auto coinciding_vertex :
//            coinciding_vertex_groups.at(coinciding_vertex_key))
//         {
//           // Skip the current vertex
//           if (coinciding_vertex == vertex_id)
//             continue;
//
//           for (const auto &neighbor_cell : v_to_c[coinciding_vertex])
//             {
//               unique_neighbors.insert(neighbor_cell);
//             }
//         }
//     }
//
//   periodic_neighbor_list.assign(unique_neighbors.begin(),
//                                 unique_neighbors.end());
// }

template void
find_cell_neighbors<2, false>(
  const parallel::distributed::Triangulation<2> &triangulation,
  dem_data_structures<2>::cells_neighbor_list   &cells_local_neighbor_list,
  dem_data_structures<2>::cells_neighbor_list   &cells_ghost_neighbor_list);

template void
find_cell_neighbors<3, false>(
  const parallel::distributed::Triangulation<3> &triangulation,
  dem_data_structures<3>::cells_neighbor_list   &cells_local_neighbor_list,
  dem_data_structures<3>::cells_neighbor_list   &cells_ghost_neighbor_list);

template void
find_cell_neighbors<2, true>(
  const parallel::distributed::Triangulation<2> &triangulation,
  dem_data_structures<2>::cells_neighbor_list   &cells_local_neighbor_list,
  dem_data_structures<2>::cells_neighbor_list   &cells_ghost_neighbor_list);

template void
find_cell_neighbors<3, true>(
  const parallel::distributed::Triangulation<3> &triangulation,
  dem_data_structures<3>::cells_neighbor_list   &cells_local_neighbor_list,
  dem_data_structures<3>::cells_neighbor_list   &cells_ghost_neighbor_list);

template void
find_cell_periodic_neighbors<2>(
  const parallel::distributed::Triangulation<2> &triangulation,
  const typename dem_data_structures<2>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  const typename dem_data_structures<2>::cell_touch_boundary_id
                                         &cell_to_pbc_mesh_id_set,
  const PeriodicBoundariesManipulator<2> &periodic_boundaries_object,
  std::vector<typename dem_data_structures<2>::cells_neighbor_list>
    &cells_local_periodic_neighbor_lists,
  std::vector<typename dem_data_structures<2>::cells_neighbor_list>
    &cells_ghost_periodic_neighbor_lists,
  std::vector<typename dem_data_structures<2>::cells_neighbor_list>
    &cells_ghost_local_periodic_neighbor_lists);

template void
find_cell_periodic_neighbors<3>(
  const parallel::distributed::Triangulation<3> &triangulation,
  const typename dem_data_structures<3>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  const typename dem_data_structures<3>::cell_touch_boundary_id
                                         &cell_to_pbc_mesh_id_set,
  const PeriodicBoundariesManipulator<3> &periodic_boundaries_object,
  std::vector<typename dem_data_structures<3>::cells_neighbor_list>
    &cells_local_periodic_neighbor_lists,
  std::vector<typename dem_data_structures<3>::cells_neighbor_list>
    &cells_ghost_periodic_neighbor_lists,
  std::vector<typename dem_data_structures<3>::cells_neighbor_list>
    &cells_ghost_local_periodic_neighbor_lists);

template void
find_full_cell_neighbors<2>(
  const parallel::distributed::Triangulation<2>     &triangulation,
  dem_data_structures<2>::cells_total_neighbor_list &cells_total_neighbor_list);

template void
find_full_cell_neighbors<3>(
  const parallel::distributed::Triangulation<3>     &triangulation,
  dem_data_structures<3>::cells_total_neighbor_list &cells_total_neighbor_list);
