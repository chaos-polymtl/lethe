#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/base/point.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/particle_insertion.h"
#include "dem/particle_wall_contact_detection.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim, dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int                cellNum = tr.n_active_cells();
  MappingQ<dim, dim> mapping(1);

  const unsigned int n_properties = 24;
  float x_min = -0.05; float y_min = -0.05; float z_min = -0.05;
  float x_max = 0.05;  float y_max = 0.05;  float z_max = 0.05;
  double dp = 0.005;
  int nInsert = 10;
  int rhop = 2500;
  Point<dim> g = {0, 0, -9.81};

  Particles::ParticleHandler<dim, dim> particle_handler(tr,
                                                        mapping,
                                                        n_properties);

  int                     nPart = 0;
  Particles::PropertyPool pool(n_properties);
  Particles::Particle<3>  particle;
   ParticleInsertion<dim, dim> ins1(x_min, y_min, z_min, x_max, y_max, z_max, dp, nInsert);
   ins1.uniformInsertion(
     particle_handler, tr, nPart,  pool, x_min, y_min, z_min, x_max, y_max, z_max, dp, nInsert, rhop, g);

  int i = 1;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++i)
    {
      deallog << "Particle " << i
              << " is inserted at: " << particle->get_location()[0] << " "
              << particle->get_location()[1] << " "
              << particle->get_location()[2] << " " << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
