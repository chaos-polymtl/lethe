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

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef FINDCELLNEIGHBORS_H_
#  define FINDCELLNEIGHBORS_H_

/**
 * Finds the neighbors lists of all the active cells in the input triangulation.
 * It is written to avoid any repetition, for instance if cell B is recognized
 * as the neighbor of cell A once, cell A will not appear in the neighbor list
 * of cell B again.
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
   * Find the neighbor list of all the active cells in the triangulation
   *
   * @param triangulation Triangulation to access the information of the cells
   * @return A vector (with size of the cell number) of sets (adjacent cells
   * of each cell). First element of each set shows the main cell itself
   */

  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
  find_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation);
};

#endif /* FINDCELLNEIGHBORS_H_ */
