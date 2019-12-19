#include <deal.II/base/point.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/contact_search.h"
#include "dem/particle_insertion.h"

using namespace dealii;

template <int dim>
void
test()
{
  Triangulation<dim, dim>              tr;
  Particles::ParticleHandler<dim, dim> particle_handler;
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

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means they are inserted in three adjacent cells in x direction
  Point<3> position1 = {-0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.4, 0, 0};
  int      id2       = 1;
  Point<3> position3 = {0.8, 0, 0};
  int      id3       = 2;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim, dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim, dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim, dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());
  Particles::ParticleIterator<dim, dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);

  Particles::Particle<dim> particle3(position3, position3, id3);
  typename Triangulation<dim, dim>::active_cell_iterator cell3 =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());
  Particles::ParticleIterator<dim, dim> pit3 =
    particle_handler.insert_particle(particle3, cell3);

  std::vector<std::pair<Particles::ParticleIterator<dim, dim>,
                        Particles::ParticleIterator<dim, dim>>>
    pairs;
  pairs = cs1.findContactPairs(particle_handler,
                               tr,
                               cellNeighbor.second,
                               cellNeighbor.first);

  for (unsigned int i = 0; i != pairs.size(); i++)
    {
      deallog << "A pair is detected: particle " << pairs[i].first->get_id()
              << " and particle " << pairs[i].second->get_id() << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  test<3>();
}
