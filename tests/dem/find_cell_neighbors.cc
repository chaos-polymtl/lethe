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

// Generate a triangulation that is twiced refined and check if the cells
// neighbors are the correct ones

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/contact_search.h"

using namespace dealii;

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, dim> tr;
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int cellNum = tr.n_active_cells();
  std::pair<std::vector<
              std::set<typename Triangulation<dim, dim>::active_cell_iterator>>,
            std::vector<typename Triangulation<dim, dim>::active_cell_iterator>>
    cellNeighbor;

  ContactSearch<dim, spacedim> cs1;
  cellNeighbor = cs1.findCellNeighbors(tr);

  int i = 0;
  for (auto cell = tr.begin_active(); cell != tr.end(); ++cell)
    {
      deallog << "neighbors of cell " << cell << " are: ";
      for (auto it = (cellNeighbor.first)[i].begin();
           it != (cellNeighbor.first)[i].end();
           ++it)
        {
          deallog << " " << *it;
        }


      deallog << std::endl;

      i++;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  test<3, 3>();
}
