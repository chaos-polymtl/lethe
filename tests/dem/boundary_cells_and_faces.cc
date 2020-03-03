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

// Check if a particle is indeed on the boundary of the system

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
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  int numRef = 2;
  triangulation.refine_global(numRef);

  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;

  FindBoundaryCellsInformation<dim> boundary_cell_object;

  boundary_cells_information =
    boundary_cell_object.find_boundary_cells_information(triangulation);

  int i = 0;
  for (unsigned int i = 0; i != boundary_cells_information.size(); ++i)
    {
      deallog << "Cell " << (boundary_cells_information[i]).cell
              << " is on system boundaries (boundary"
              << (boundary_cells_information[i]).boundary_id << ")"
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
