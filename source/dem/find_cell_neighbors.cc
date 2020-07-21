#include <dem/find_cell_neighbors.h>

using namespace dealii;

// The constructor of this class is empty
template <int dim> FindCellNeighbors<dim>::FindCellNeighbors() {}

// This function finds the neighbor list of all the active cells in the
// triangulation
template <int dim>
void FindCellNeighbors<dim>::find_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
        &cells_local_neighbor_list,
    std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
        &cells_ghost_neighbor_list) {
  // Number of active cells in the triangulation
  // int cell_number = triangulation.n_active_cells();

  // The output vector of the function. It is a vector with size of the number
  // of active cells. Each element of this vector is a set which shows the
  // corresponding adjacent cells of the main cell. The first element of the set
  // is the main cell.map

  std::vector<typename Triangulation<dim>::active_cell_iterator>
      local_neighbor_vector;
  std::vector<typename Triangulation<dim>::active_cell_iterator>
      ghost_neighbor_vector;
  // This vector is used to avoid repetition of adjacent cells. For instance if
  // cell B is recognized as the neighbor of cell A, cell A will not be added to
  // the neighbor list of cell B again. This is done using the totall_cell_list
  // vector
  std::vector<typename Triangulation<dim>::active_cell_iterator>
      totall_cell_list;

  // Cell iterator counter
  // int cell_number_iterator = 0;

  // For each cell, the cell vertices are found and used to find adjacent cells.
  // The reason is to find the cells located on the corners of the main cell.
  auto v_to_c = GridTools::vertex_to_cell_map(triangulation);

  unsigned int cell_number_iterator = 0;
  // Looping over cells
  for (const auto &cell : triangulation.active_cell_iterators()) {

    if (cell->is_locally_owned()) {

      // The first element of each set (eacmpirun -np 2 ./dem_2d
      // packing_in_square.prm h element of the vector) is the cell itself.
      local_neighbor_vector.push_back(cell);

      totall_cell_list.push_back(cell);

      for (unsigned int vertex = 0;
           vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex) {
        for (const auto &neighbor : v_to_c[cell->vertex_index(vertex)]) {

          if (neighbor->is_locally_owned()) {
            auto search_iterator = std::find(totall_cell_list.begin(),
                                             totall_cell_list.end(), neighbor);
            auto local_search_iterator =
                std::find(local_neighbor_vector.begin(),
                          local_neighbor_vector.end(), neighbor);

            // If the cell (neighbor) is not present in the
            // total_cell_list vector, it will be added as the neighbor of
            // the main cell
            // ("cell") and also to the total_cell_list to avoid
            // repetition for next cells.
            if (search_iterator == totall_cell_list.end() &&
                local_search_iterator == local_neighbor_vector.end()) {
              local_neighbor_vector.push_back(neighbor);
            }

          } else if (neighbor->is_ghost()) {
            auto ghost_search_iterator =
                std::find(ghost_neighbor_vector.begin(),
                          ghost_neighbor_vector.end(), neighbor);
            if (ghost_search_iterator == ghost_neighbor_vector.end()) {
              if (ghost_neighbor_vector.empty()) {
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
    ++cell_number_iterator;
  }
}

template class FindCellNeighbors<2>;
template class FindCellNeighbors<3>;
