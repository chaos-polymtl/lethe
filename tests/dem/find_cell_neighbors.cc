#include "dem/contact_search.h"
#include "../tests.h"
#include <iostream>
#include <vector>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>

using namespace dealii;

int main(int argc, char **argv)
{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(
      argc, argv, numbers::invalid_unsigned_int);
    parallel::distributed::Triangulation<3, 3> tr(MPI_COMM_WORLD);
    GridGenerator::hyper_cube(tr, -1, 1, true);
    int numRef = 2;
    tr.refine_global(numRef);
    initlog();
 int cellNum = tr.n_active_cells();
    std::pair <std::vector<std::set<Triangulation<3>::active_cell_iterator>>,std::vector<Triangulation<3>::active_cell_iterator>> cellNeighbor;

  ContactSearch cs1;
  cellNeighbor=cs1.findCellNeighbors(cellNum, tr);

int i = 0;
  for (Triangulation<3>::active_cell_iterator cell = tr.begin_active();
       cell != tr.end();
       ++cell)
{
deallog << "neighbors of cell " << cell << " are: ";
            for (auto it=(cellNeighbor.first)[i].begin(); it != (cellNeighbor.first)[i].end(); ++it)
{
     deallog << " " << *it;
}


deallog << std::endl;

        i++;
}

}
