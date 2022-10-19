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

// This function finds the full neighbor list (with repetition) of all the
// active cells in the triangulation. We need this function for the particle-
// floating mesh contact force calculations. In particle-floating mesh contacts
// , we need all the particles located in ALL the neighbor cells of the main
// cell to search for possible collisions with the floating mesh.
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
