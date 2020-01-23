#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/contact_search.h"
#include "dem/integration.h"
#include "dem/particle_insertion.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim, dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  MappingQ<dim, dim> mapping(1);

  const unsigned int n_properties = 24;
  Point<dim>         g            = {0, 0, -9.81};
  float              dt           = 0.00001;

  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);



  // inserting one particle at x = 0 , y = 0 and z = 0 m
  // initial velocity of particles = 0, 0, 0 m/s
  // gravitational acceleration = 0, 0, -9.81 m/s2
  Point<3> position1 = {0, 0, 0};
  int      id1       = 0;

  Particles::Particle<dim> particle1(position1, position1, id1);


  typename Triangulation<dim, dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());


  Particles::ParticleIterator<dim, dim> pit =
    particle_handler.insert_particle(particle1, cell1);



  pit->get_properties()[0] = id1;

  pit->get_properties()[1] = 1;
  pit->get_properties()[2] = 0.005;
  pit->get_properties()[3] = 2500;
  // Position
  pit->get_properties()[4] = 0;
  pit->get_properties()[5] = 0;
  pit->get_properties()[6] = 0;
  // Velocity
  pit->get_properties()[7] = 0;
  pit->get_properties()[8] = 0;
  pit->get_properties()[9] = 0;
  // Acceleration
  pit->get_properties()[10] = 0;
  pit->get_properties()[11] = 0;
  pit->get_properties()[12] = -9.81;
  // Force
  pit->get_properties()[13] = 0;
  pit->get_properties()[14] = 0;
  pit->get_properties()[15] = 0;
  // w
  pit->get_properties()[16] = 0;
  pit->get_properties()[17] = 0;
  pit->get_properties()[18] = 0;
  // mass and moment of inertia
  pit->get_properties()[19] = 1;
  pit->get_properties()[20] = 1;


  Integration<dim> Integ1;
  Integ1.velocity_verlet_integration(particle_handler, g, dt);

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      deallog << "The new position of the particle in z direction after " << dt
              << " seconds is: " << particle->get_properties()[6] << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
