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
#include "dem/particle_wall_contact_detection.h"


using namespace dealii;

template <int dim>
void
test()
{
  Triangulation<dim, dim> tr;
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  std::vector<std::tuple<int,
                         Triangulation<3>::active_cell_iterator,
                         int,
                         Point<3>,
                         Point<3>>>
    boundaryCellInfo;

  ParticleWallContactDetection<dim> pw1;
  pw1.boundaryCellsAndFaces(tr, boundaryCellInfo);

  int i = 0;
  for (unsigned int i = 0; i != boundaryCellInfo.size(); ++i)
    {
      deallog << "Cell " << std::get<1>(boundaryCellInfo[i])
              << " is on system boundaries (boundary"
              << std::get<0>(boundaryCellInfo[i]) << ")" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  test<3>();
}
