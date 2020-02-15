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
#include "dem/insertion_info_struct.h"
#include "dem/physical_info_struct.h"
#include "dem/uniform_insertion.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim, dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int                           cellNum = tr.n_active_cells();
  MappingQ<dim, dim>            mapping(1);
  insertion_info_struct<dim, dim> insertion_info_struct;
  physical_info_struct<dim> physical_info_struct;

  const unsigned int n_properties               = 24;
  insertion_info_struct.x_min                   = -0.05;
  insertion_info_struct.y_min                   = -0.05;
  insertion_info_struct.z_min                   = -0.05;
  insertion_info_struct.x_max                   = 0.05;
  insertion_info_struct.y_max                   = 0.05;
  insertion_info_struct.z_max                   = 0.05;
  insertion_info_struct.inserted_number_at_step = 10;
  insertion_info_struct.distance_threshold      = 2;
  physical_info_struct.particle_diameter        = 0.005;
  physical_info_struct.particle_density         = 2500;


  Particles::ParticleHandler<dim, dim> particle_handler(tr,
                                                        mapping,
                                                        n_properties);

  Particles::PropertyPool property_pool(n_properties);
  Particles::Particle<3>  particle;

   UniformInsertion<dim, dim> ins1(physical_info_struct,
                                          insertion_info_struct);
  ins1.insert(particle_handler, tr, property_pool, physical_info_struct,
              insertion_info_struct);

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
