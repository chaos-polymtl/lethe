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

  ParticleWallContactDetection pw1;
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
