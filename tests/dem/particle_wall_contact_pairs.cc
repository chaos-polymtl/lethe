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
#include "dem/find_boundary_cells_information.h"
#include "dem/pw_broad_search.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  Particles::ParticleHandler<dim> particle_handler;
  int                             num_particles = 1;
  std::pair<
    std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>,
    std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cellNeighbor;

  // inserting three particles at x = -0.4 , x = 0.4 and x = 0.8
  // which means only particle 3 is located in a boundary cell
  Point<dim> position1 = {-0.4, 0, 0};
  int        id1       = 0;
  Point<dim> position2 = {0.4, 0, 0};
  int        id2       = 1;
  Point<dim> position3 = {0.8, 0, 0};
  int        id3       = 2;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);

  Particles::Particle<dim> particle3(position3, position3, id3);
  typename Triangulation<dim>::active_cell_iterator cell3 =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, cell3);

  std::vector<boundary_cells_info_struct<dim>> boundary_cells_information;
  FindBoundaryCellsInformation<dim>            boundary_cell_object;
  boundary_cells_information =
    boundary_cell_object.find_boundary_cells_information(tr);

  PWBroadSearch<dim> pw1;
  std::vector<
    std::tuple<std::pair<typename Particles::ParticleIterator<dim>, int>,
               Tensor<1, dim>,
               Point<dim>>>
    pwContactList(num_particles);
  pw1.find_PW_Contact_Pairs(boundary_cells_information,
                            particle_handler,
                            pwContactList);

  for (unsigned int i = 0; i != pwContactList.size(); i++)
    {
      deallog << "Particle " << std::get<0>(pwContactList[i]).first->get_id()
              << " is located in a boundary cell" << std::endl;
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
