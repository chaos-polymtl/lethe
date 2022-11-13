#include <dem/find_cell_neighbors.h>

using namespace dealii;

// The constructor of this class is empty
template <int dim>
FindCellNeighbors<dim>::FindCellNeighbors()
{}

// This function finds the neighbor list (without repetition) of all the active
// cells in the triangulation
template <int dim>
void
FindCellNeighbors<dim>::find_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list)
{
  // The output vectors of the function are cells_local_neighbor_list and
  // cells_ghost_neighbor_list. The first contains all the local neighbors cells
  // of all local cells; while the second contains all the ghost cells of all
  // local cells. They are two vectors with size of the number of active cells.
  // The first elements of all vectors are the main cells
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    local_neighbor_vector;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    ghost_neighbor_vector;

  // This vector is used to avoid repetition of adjacent cells. For instance if
  // cell B is recognized as the neighbor of cell A, cell A will not be added to
  // the neighbor list of cell B again. This is done using the total_cell_list
  // vector
  std::vector<typename Triangulation<dim>::active_cell_iterator>
    total_cell_list;

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

          total_cell_list.push_back(cell);

          for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
            {
              for (const auto &neighbor : v_to_c[cell->vertex_index(vertex)])
                {
                  if (neighbor->is_locally_owned())
                    {
                      auto search_iterator = std::find(total_cell_list.begin(),
                                                       total_cell_list.end(),
                                                       neighbor);
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
FindCellNeighbors<dim>::find_cell_periodic_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const typename DEM::dem_data_structures<dim>::cell_container
    &periodic_cells_container,
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list)
{
  typename DEM::dem_data_structures<dim>::cell_container
    local_periodic_neighbor_vector;
  typename DEM::dem_data_structures<dim>::cell_container
    ghost_periodic_neighbor_vector;

  typename DEM::dem_data_structures<dim>::cell_container total_cell_list;

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

  // Looping over cells at outlet (one side of the periodic boundaries)
  for (const auto &cell : periodic_cells_container)
    {
      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
          // The first element of each vector is the cell itself.
          local_periodic_neighbor_vector.push_back(cell);
          total_cell_list.push_back(cell);

          // Empty list of periodic cell neighbor
          typename DEM::dem_data_structures<dim>::cell_container
            periodic_neighbor_list;

          // Loop over vertices of the cell
          for (unsigned int vertex = 0; vertex < cell->n_vertices(); ++vertex)
            {
              // Get global id of vertex
              unsigned int vertex_id = cell->vertex_index(vertex);

              // Check if vertex is at periodic boundary, should be a key if so
              if (vertex_to_coinciding_vertex_group.find(vertex_id) !=
                  vertex_to_coinciding_vertex_group.end())
                {
                  // Get the coinciding vertex key to the current vertex
                  unsigned int coinciding_vertex_key =
                    vertex_to_coinciding_vertex_group[vertex_id];

                  // Store the neighbor cells in list
                  for (auto coinciding_vertex :
                       coinciding_vertex_groups[coinciding_vertex_key])
                    {
                      // Skip the current vertex since we want only cells linked
                      // to the periodic vertices
                      if (coinciding_vertex != vertex_id)
                        {
                          // Loop over all periodic neighbor
                          for (auto neighbor_id : v_to_c[coinciding_vertex])
                            {
                              periodic_neighbor_list.push_back(neighbor_id);
                            }
                        }
                    }
                }
            }

          for (const auto &periodic_neighbor : periodic_neighbor_list)
            {
              if (periodic_neighbor->is_locally_owned())
                {
                  auto search_iterator = std::find(total_cell_list.begin(),
                                                   total_cell_list.end(),
                                                   periodic_neighbor);
                  auto local_search_iterator =
                    std::find(local_periodic_neighbor_vector.begin(),
                              local_periodic_neighbor_vector.end(),
                              periodic_neighbor);

                  // If the cell (neighbor) is a local cell and not present
                  // in the total_cell_list vector, it will be added as the
                  // neighbor of the main cell
                  // ("cell") and also to the total_cell_list to avoid
                  // repetition for next cells.
                  if (search_iterator == total_cell_list.end() &&
                      local_search_iterator ==
                        local_periodic_neighbor_vector.end())
                    {
                      local_periodic_neighbor_vector.push_back(
                        periodic_neighbor);
                    }

                  // If the neighbor cell is a ghost, it should be added in
                  // the ghost_neighbor_vector container
                }
              else if (periodic_neighbor->is_ghost())
                {
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
    }
}

// This function finds the full neighbor list (with repetition) of all the
// active cells in the triangulation. We need this function for the particle-
// floating mesh contact force calculations. In particle-floating mesh
// contacts , we need all the particles located in ALL the neighbor cells of
// the main cell to search for possible collisions with the floating mesh.
template <int dim>
void
FindCellNeighbors<dim>::find_full_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
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

template class FindCellNeighbors<2>;
template class FindCellNeighbors<3>;
