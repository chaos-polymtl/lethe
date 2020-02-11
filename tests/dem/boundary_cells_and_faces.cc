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

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/find_boundary_cells_information.h"


using namespace dealii;

template <int dim>
void
test()
{
  Triangulation<dim, dim> tr;
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  std::vector<boundary_cells_info_struct<dim>> boundaryCellInfo;


 FindBoundaryCellsInformation<dim, dim> pw1;
  boundaryCellInfo = pw1.find_boundary_cells_information(tr);

  int i = 0;
  for (unsigned int i = 0; i != boundaryCellInfo.size(); ++i)
    {
      deallog << "Cell " << (boundaryCellInfo[i]).cell
              << " is on system boundaries (boundary"
              << (boundaryCellInfo[i]).boundary_id << ")" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  test<3>();
}
