// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
  // local cells. They are two vectors with size of the number of active cells.
  // The first elements of all vectors are the main cells
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
  const typename dem_data_structures<dim>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  typename dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  typename dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list)
{
  typename dem_data_structures<dim>::cell_vector local_periodic_neighbor_vector;
  typename dem_data_structures<dim>::cell_vector ghost_periodic_neighbor_vector;
  typename dem_data_structures<dim>::cell_vector
    ghost_local_periodic_neighbor_vector;

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

  // Looping over cells struct at periodic boundaries 0
  for (const auto &pb_cell_struct : periodic_boundaries_cells_information)
    {
      auto &cell = pb_cell_struct.second.cell;

      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
          // The first element of each vector is the cell itself.
          local_periodic_neighbor_vector.push_back(cell);

          // Store every main cell
          total_cell_list.insert(cell);

          // Empty list of periodic cell neighbor
          typename dem_data_structures<dim>::cell_vector periodic_neighbor_list;

          // Get the periodic neighbor of the cell
          get_periodic_neighbor_list<dim>(cell,
                                          coinciding_vertex_groups,
                                          vertex_to_coinciding_vertex_group,
                                          v_to_c,
                                          periodic_neighbor_list);

          for (const auto &periodic_neighbor : periodic_neighbor_list)
            {
              if (periodic_neighbor->is_locally_owned())
                {
                  // Check if the neighbor cell has been processed as a
                  // main cell
                  auto search_iterator =
                    total_cell_list.find(periodic_neighbor);

                  // Check if the neighbor cell is in already in the
                  // local_periodic_neighbor_vector
                  // Note from Gabo: I don't understand how this could happen
                  // since we are looping over cell on boundary 1.
                  // I think the only case would be if periodic_neighbor_list
                  // has duplicate cell. If this is the case, it is weird that
                  // get_periodic_neighbor_list returns a list with duplicates.
                  auto local_search_iterator =
                    std::find(local_periodic_neighbor_vector.begin(),
                              local_periodic_neighbor_vector.end(),
                              periodic_neighbor);

                  // If the cell neighbor is a local cell and not present
                  // in the total_cell_list vector, it will be added as the
                  // neighbor of the main cell.
                  if (search_iterator == total_cell_list.end() &&
                      local_search_iterator ==
                        local_periodic_neighbor_vector.end())
                    {
                      local_periodic_neighbor_vector.push_back(
                        periodic_neighbor);
                    }
                }
              else if (periodic_neighbor->is_ghost())
                {
                  // If the cell neighbor is a ghost, it should be added in
                  // the ghost_periodic_neighbor_vector container
                  auto ghost_search_iterator =
                    std::find(ghost_periodic_neighbor_vector.begin(),
                              ghost_periodic_neighbor_vector.end(),
                              periodic_neighbor);
                  if (ghost_search_iterator ==
                      ghost_periodic_neighbor_vector.end())
                    {
                      if (ghost_periodic_neighbor_vector.empty())
                        {
                          ghost_periodic_neighbor_vector.push_back(cell);
                        }

                      ghost_periodic_neighbor_vector.push_back(
                        periodic_neighbor);
                    }
                }
            }
          if (!local_periodic_neighbor_vector.empty())
            cells_local_periodic_neighbor_list.push_back(
              local_periodic_neighbor_vector);
          if (!ghost_periodic_neighbor_vector.empty())
            cells_ghost_periodic_neighbor_list.push_back(
              ghost_periodic_neighbor_vector);
          local_periodic_neighbor_vector.clear();
          ghost_periodic_neighbor_vector.clear();
        }
      else if (cell->is_ghost())
        {
          // Since periodic cells are mapped on one side only (cells on pb 0
          // with cells on pb 1), we need a 3rd container for ghost-local
          // contacts for force calculation. Here we store local neighbors of
          // ghost cells.

          // The first element of each vector is the cell itself
          ghost_local_periodic_neighbor_vector.push_back(cell);
          total_ghost_cell_list.insert(cell);

          // Empty list of periodic cell neighbor
          typename dem_data_structures<dim>::cell_vector periodic_neighbor_list;

          // Get the periodic neighbor of the cell
          get_periodic_neighbor_list<dim>(cell,
                                          coinciding_vertex_groups,
                                          vertex_to_coinciding_vertex_group,
                                          v_to_c,
                                          periodic_neighbor_list);

          for (const auto &periodic_neighbor : periodic_neighbor_list)
            {
              if (periodic_neighbor->is_locally_owned())
                {
                  auto search_iterator =
                    total_ghost_cell_list.find(periodic_neighbor);

                  auto local_search_iterator =
                    std::find(ghost_local_periodic_neighbor_vector.begin(),
                              ghost_local_periodic_neighbor_vector.end(),
                              periodic_neighbor);

                  if (search_iterator == total_ghost_cell_list.end() &&
                      local_search_iterator ==
                        ghost_local_periodic_neighbor_vector.end())
                    {
                      ghost_local_periodic_neighbor_vector.push_back(
                        periodic_neighbor);
                    }
                }
            }
          if (!ghost_local_periodic_neighbor_vector.empty())
            cells_ghost_local_periodic_neighbor_list.push_back(
              ghost_local_periodic_neighbor_vector);
          ghost_local_periodic_neighbor_vector.clear();
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

template <int dim>
void
get_periodic_neighbor_list(
  const typename Triangulation<dim>::active_cell_iterator &cell,
  const std::map<unsigned int, std::vector<unsigned int>>
                                             &coinciding_vertex_groups,
  const std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group,
  const std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
                                                 &v_to_c,
  typename dem_data_structures<dim>::cell_vector &periodic_neighbor_list)
{
  // Loop over all vertices of the cell (periodic and non periodic)
  for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
    {
      // Get global id of vertex
      unsigned int vertex_id = cell->vertex_index(vertex);

      // Check if vertex is at periodic boundary, there should be a key if so
      if (vertex_to_coinciding_vertex_group.find(vertex_id) !=
          vertex_to_coinciding_vertex_group.end())
        {
          // Get the coinciding vertex key to the current vertex
          unsigned int coinciding_vertex_key =
            vertex_to_coinciding_vertex_group.at(vertex_id);

          // Store the neighbor cells in list
          for (auto coinciding_vertex :
               coinciding_vertex_groups.at(coinciding_vertex_key))
            {
              // Skip the current vertex since we want only cells linked
              // to the periodic vertices
              if (coinciding_vertex != vertex_id)
                {
                  // Loop over all periodic neighbor
                  for (const auto &neighbor_id : v_to_c[coinciding_vertex])
                    {
                      periodic_neighbor_list.push_back(neighbor_id);
                    }
                }
            }
        }
    }
}

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
  const dem_data_structures<2>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  dem_data_structures<2>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  dem_data_structures<2>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list);

template void
find_cell_periodic_neighbors<3>(
  const parallel::distributed::Triangulation<3> &triangulation,
  const dem_data_structures<3>::periodic_boundaries_cells_info
    &periodic_boundaries_cells_information,
  dem_data_structures<3>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  dem_data_structures<3>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list);

template void
find_full_cell_neighbors<2>(
  const parallel::distributed::Triangulation<2>     &triangulation,
  dem_data_structures<2>::cells_total_neighbor_list &cells_total_neighbor_list);

template void
find_full_cell_neighbors<3>(
  const parallel::distributed::Triangulation<3>     &triangulation,
  dem_data_structures<3>::cells_total_neighbor_list &cells_total_neighbor_list);
