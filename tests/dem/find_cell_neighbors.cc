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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

// This test generates a triangulation that is twiced refined and check if the
// cells neighbors are the correct ones

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <dem/find_cell_neighbors.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim> void test() {
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation, -1 * hyper_cube_length,
                            hyper_cube_length, true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);

  // Finding the cell neighbors
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_local_neighbor_list;
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
      cells_ghost_neighbor_list;

  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_object.find_cell_neighbors(
      triangulation, cells_local_neighbor_list, cells_ghost_neighbor_list);

  // Output
  int i = 0;
  for (auto cell = triangulation.begin_active(); cell != triangulation.end();
       ++cell) {
    deallog << "neighbors of cell " << cell << " are: ";
    for (auto iterator = cells_local_neighbor_list[i].begin();
         iterator != cells_local_neighbor_list[i].end(); ++iterator) {
      deallog << " " << *iterator;
    }
    deallog << std::endl;

    ++i;
  }
}

int main(int argc, char **argv) {
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
      argc, argv, numbers::invalid_unsigned_int);
  test<3>();
}
