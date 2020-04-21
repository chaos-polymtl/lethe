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

// This test finds the boundary cells in a system and reports the corresponding
// information of these cells

#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <dem/find_boundary_cells_information.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the triangulation and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);

  // Fining boundary cellds information
  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;
  FindBoundaryCellsInformation<dim>            boundary_cell_object;
  boundary_cells_information =
    boundary_cell_object.find_boundary_cells_information(triangulation);

  // Reporting the information of boundary cells
  for (auto boundary_cells_information_iterator =
         boundary_cells_information.begin();
       boundary_cells_information_iterator != boundary_cells_information.end();
       ++boundary_cells_information_iterator)
    {
      deallog << "Cell " << boundary_cells_information_iterator->cell
              << " is on system boundaries (boundary"
              << boundary_cells_information_iterator->boundary_id << ")"
              << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);
  test<3>();
}
