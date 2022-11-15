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
   * Finds the neighbor list (without repetition) of all the active cells in the
   * triangulation
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
   * Finds the neighbor list (without repetition) of all the active cells in the
   * triangulation
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param cells_local_periodic_neighbor_list A vector (with size of the local cell
   * number) of vectors (local adjacent cells of each local cell). First element
   * of each set shows the main cell itself
   * @param cells_ghost_periodic_neighbor_list A vector (with size of the local cell
   * number) of vectors (ghost adjacent cells of each local cell). First element
   * of each set shows the main cell itself
   */

  void
  find_cell_periodic_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_local_periodic_neighbor_list,
    typename DEM::dem_data_structures<dim>::cells_neighbor_list
      &cells_ghost_periodic_neighbor_list);

  /**
   * Finds the full neighbor list (with repetition) of all the active cells in
   * the triangulation. This function is used in particle-floating mesh
   * contacts, as in this situation, we need to define all the particles located
   * in the neighbor cells of the background cell (cut by the floating mesh) as
   * contact candidates
   *
   * @param triangulation Triangulation to access the information of the cells
   * @param cells_total_neighbor_list An unordered_map (with size of the local cell
   * number) of vectors (all adjacent cells of each local cell)
   */

  void
  find_full_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
      &cells_total_neighbor_list);
};

#endif /* find_cell_neighbors_h */
