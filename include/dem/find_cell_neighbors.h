/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <dem/data_containers.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>


using namespace dealii;

#ifndef find_cell_neighbors_h
#  define find_cell_neighbors_h

/**
 * Finds the neighbors lists of all the active cells in the input triangulation.
 * find_cell_neighbors() is written to avoid any repetition, for instance if
 * cell B is recognized as the neighbor of cell A once, cell A will not appear
 * in the neighbor list of cell B again. On the other hand,
 * find_full_cell_neighbors() function finds the neighbors list with this
 * repetition.
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class FindCellNeighbors
{
public:
  FindCellNeighbors<dim>();

  /**
   * @brief Finds the neighbor list (without repetition) of all the active
   * cells in the triangulation. It gets the vertices of the cells to get lists
   * of the neighbor for each cells. There is some check to prevent repetition
   * of a cell in a list (up to 4 vertices can have the same cell in common in
   * 3D). 2 types of container are used for cell neighbors : local-local cells
   * and local-ghost cells.
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param cells_local_neighbor_list A vector (with size of the local cell
   * number) of vectors (local adjacent cells of each local cell). First element
   * of each set shows the main cell itself
   * @param cells_ghost_neighbor_list A vector (with size of the local cell
   * number) of vectors (ghost adjacent cells of each local cell). First element
   * of each set shows the main cell itself
   */

  void
  find_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_local_neighbor_list,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_ghost_neighbor_list);

  /**
   * @brief Finds the periodic neighbor list (without repetition) of all the
   * active cells in the triangulation. It gets the coinciding vertices of the
   * cells at the periodic boundary 0 to get lists of the periodic neighbor
   * cells on the periodic boundary 1 for each cells. There is some check to
   * prevent repetition of a cell in a list (up to 4 vertices can have the
   * same cell in common in 3D). Also, 3 types of container are used for
   * periodic mapping of cell neighbors : local-local cells, local-ghost cells
   * and ghost-local cells. The last container is necessary since the mapping
   * are only from periodic boundary 0 to periodic boundary 1 and the
   * ghost-local particle pairs need distinction for proper handling of search
   * of particle pairs and contact forces.
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param periodic_boundaries_cells_information A container of information
   * related to the pairs of cell at periodic boundaries, used to get the cells
   * on periodic boundary 0
   * @param cells_local_periodic_neighbor_list A vector (with size of the local
   * cell number) of vectors (local adjacent cells of each local cell). First
   * element of each set shows the main cell itself
   * @param cells_ghost_periodic_neighbor_list A vector (with size of the local
   * cell number) of vectors (ghost adjacent cells of each local cell). First
   * element of each set shows the main cell itself
   * @param cells_ghost_local_periodic_neighbor_list A vector (with size of the
   * ghost cell number) of vectors (local adjacent cells of each ghost cell).
   * First element of each set shows the main ghost cell itself
   */

  void
  find_cell_periodic_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_local_periodic_neighbor_list,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_ghost_periodic_neighbor_list,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_ghost_local_periodic_neighbor_list);

  /**
   * Finds the full neighbor list (with repetition) of all the active cells in
   * the triangulation. This function is used in particle-floating mesh
   * contacts, as in this situation, we need to define all the particles located
   * in the neighbor cells of the background cell (cut by the floating mesh) as
   * contact candidates
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param cells_total_neighbor_list An unordered_map (with size of the local
   * cell number) of vectors (all adjacent cells of each local cell)
   */

  void
  find_full_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
      &cells_total_neighbor_list);

private:
  /**
   * @brief Generate a periodic neighbor cells list of the cell. With the
   * coinciding_vertex_groups and the vertex_to_coinciding_vertex_group, we can
   * get the coinciding vertices on a periodic boundary and with v_to_c we can
   * get the list of the periodic neighbor cells to the main cell.
   *
   * Here is an example to understand how the search works :
   * Cell 8 have vertices 0, 1, 2, 3 on the periodic boundary. First,
   * it checks in the coinciding_vertex_groups if the vertex 0 is a key in the
   * map, and if it is, it gets a label, let's say label is 10. In the
   * vertex_to_coinciding_vertex_group, it can finds the the coinciding vertices
   * with the label. Using the label 10, it gets the vertices 0 and 34, which
   * means that vertices 0 and 34 are periodic. We do not want the cells
   * attached to the vertex 0 since they are already found with the regular
   * find_cell_neighbor function, so it skips vertex 0, but gets the cells 21,
   * 22, 23, 24 attached to the vertex 34. The same checks are done for vertices
   * 1, 2 and 3 and all periodic cells are return as a vector.
   *
   * @param cell The cell that needs the periodic neighbor list
   * @param coinciding_vertex_groups A map of coinciding vertices labeled by an
   * arbitrary element from them
   * @param vertex_to_coinciding_vertex_group Map of a vertex to the label of a
   * group of coinciding vertices
   * @param v_to_c A vector of set with adjacent cells of all the vertices
   * @param periodic_neighbor_list A vector which is the list of periodic cell
   * neighbors
   */

  void
  get_periodic_neighbor_list(
    const typename Triangulation<dim>::active_cell_iterator &cell,
    const std::map<unsigned int, std::vector<unsigned int>>
      &coinciding_vertex_groups,
    const std::map<unsigned int, unsigned int>
      &vertex_to_coinciding_vertex_group,
    const std::vector<
      std::set<typename Triangulation<dim>::active_cell_iterator>> &v_to_c,
    typename DEM::dem_data_structures<dim>::cell_container
      &periodic_neighbor_list);
};

#endif /* find_cell_neighbors_h */
