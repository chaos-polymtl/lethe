#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/contact_search.h"

using namespace dealii;

template <int dim>
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

  ContactSearch cs1;
  cellNeighbor = cs1.findCellNeighbors(cellNum, tr);

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
  test<3>();
}
