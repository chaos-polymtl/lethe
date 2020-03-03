#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/dem_solver_parameters.h"
#include "dem/uniform_insertion.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int                      cellNum = tr.n_active_cells();
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  const unsigned int n_properties                 = 21;
  dem_parameters.insertionInfo.x_min              = -0.05;
  dem_parameters.insertionInfo.y_min              = -0.05;
  dem_parameters.insertionInfo.z_min              = -0.05;
  dem_parameters.insertionInfo.x_max              = 0.05;
  dem_parameters.insertionInfo.y_max              = 0.05;
  dem_parameters.insertionInfo.z_max              = 0.05;
  dem_parameters.insertionInfo.inserted_this_step = 10;
  dem_parameters.insertionInfo.distance_threshold = 2;
  dem_parameters.physicalProperties.diameter      = 0.005;
  dem_parameters.physicalProperties.density       = 2500;

  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  Particles::PropertyPool property_pool(n_properties);
  Particles::Particle<3>  particle;

  UniformInsertion<dim> ins1(dem_parameters);
  ins1.insert(particle_handler, tr, property_pool, dem_parameters);

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
