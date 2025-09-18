// SPDX-FileCopyrightText: Copyright (c) 2020- 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_find_cell_neighbors_h
#define lethe_find_cell_neighbors_h

#include <dem/data_containers.h>

#include <deal.II/distributed/tria.h>

using namespace dealii;

/**
 * @brief Finds the neighbor list (without repetition) of all the active
 * cells in the triangulation. It gets the vertices of the cells to get lists
 * of the neighbor for each cell. There is some check to prevent repetition
 * of a cell in a list (up to 8 vertices can have the same cell in common in
 * 3D).
 * 2 types of container are used for cell neighbors : local-local cells
 * and local-ghost cells.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the problem is solved.
 * @tparam reciprocal A boolean that denotes if the cells_local_neighbor_list of
 * cell i and j contains cell j and i, respectively.
 *
 * @param triangulation Triangulation to access the information of the cells
 * @param cells_local_neighbor_list A vector (with size of the local cell
 * number) of vectors (local adjacent cells of each local cell). First element
 * of each set shows the main cell itself
 * @param cells_ghost_neighbor_list A vector (with size of the local cell
 * number) of vectors (ghost adjacent cells of each local cell). First element
 * of each set shows the main cell itself
 */
template <int dim, bool reciprocal = false>
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
 * cells on the periodic boundary 1 for each cell. There is some check to
 * prevent repetition of a cell in a list (up to 8 vertices can have the
 * same cell in common in 3D).
 * 3 types of container are used for periodic mapping of cell neighbors :
 * local-local cells, local-ghost cells and ghost-local cells. The last
 * container is necessary since the mapping are only from periodic boundary 0
 * to periodic boundary 1 and the ghost-local particle pairs need distinction
 * for proper handling of search of particle pairs and contact forces.
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
template <int dim>
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
 * @brief Finds the full neighbor list (with repetition) of all the active
 * cells in the triangulation. This function is used in particle-floating mesh
 * contacts, as in this situation, we need to define all the particles located
 * in the neighbor cells of the background cell (cut by the floating mesh) as
 * contact candidates
 *
 * @param triangulation Triangulation to access the information of the cells
 * @param cells_total_neighbor_list An unordered_map (with size of the local
 * cell number) of vectors (all adjacent cells of each local cell)
 */
template <int dim>
void
find_full_cell_neighbors(
  const parallel::distributed::Triangulation<dim> &triangulation,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
    &cells_total_neighbor_list);

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
 * vertex_to_coinciding_vertex_group, it can find the coinciding vertices
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
template <int dim>
void
get_periodic_neighbor_list(
  const typename Triangulation<dim>::active_cell_iterator &cell,
  const std::map<unsigned int, std::vector<unsigned int>>
                                             &coinciding_vertex_groups,
  const std::map<unsigned int, unsigned int> &vertex_to_coinciding_vertex_group,
  const std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
                                                      &v_to_c,
  typename DEM::dem_data_structures<dim>::cell_vector &periodic_neighbor_list);

#endif
