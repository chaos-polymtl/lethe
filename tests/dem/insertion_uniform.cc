#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>

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

  std::string      filename = "../dem.prm";
  ParameterHandler prm;
  ParametersDEM<3> DEMparam;
  DEMparam.declare(prm);
  prm.parse_input(filename);
  DEMparam.parse(prm);

  const unsigned int n_properties = DEMparam.outputProperties.numProperties;

  Particles::ParticleHandler<dim, dim> particle_handler(tr,
                                                        mapping,
                                                        n_properties);

  int                     nPart = 0;
  Particles::PropertyPool pool(DEMparam.outputProperties.numProperties);
  Particles::Particle<3>  particle;
  ParticleInsertion       ins1(DEMparam);
  ins1.uniformInsertion(particle_handler, tr, DEMparam, nPart, pool, particle);

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
