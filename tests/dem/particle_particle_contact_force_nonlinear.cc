#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_solver_parameters.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/find_cell_neighbors.h"
#include "dem/pp_broad_search.h"
#include "dem/pp_fine_search.h"
#include "dem/pp_nonlinear_force.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, -1, 1, true);
  int numRef = 2;
  triangulation.refine_global(numRef);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  const unsigned int n_properties = 21;
  Tensor<1, dim>     g{{0, 0, -9.81}};
  float              dt = 0.00001;

  double particle_diameter                                           = 0.005;
  int    particle_density                                            = 2500;
  dem_parameters.physicalProperties.Youngs_modulus_particle          = 50000000;
  dem_parameters.physicalProperties.Poisson_ratio_particle           = 0.3;
  dem_parameters.physicalProperties.restitution_coefficient_particle = 0.5;
  dem_parameters.physicalProperties.friction_coefficient_particle    = 0.5;
  dem_parameters.physicalProperties.rolling_friction_particle        = 0.1;

  Particles::ParticleHandler<dim> particle_handler(triangulation,
                                                   mapping,
                                                   n_properties);

  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    cellNeighbor;

  FindCellNeighbors<dim> cn1;
  cellNeighbor = cn1.find_cell_neighbors(triangulation);

  PPBroadSearch<dim> ppbs;
  PPFineSearch<dim>  ppfs;

  Point<3> position1     = {0.4, 0, 0};
  int      id1           = 0;
  Point<3> position2     = {0.40499, 0, 0};
  int      id2           = 1;
  int      num_Particles = 2;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[0]  = id1;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = particle_diameter;
  pit1->get_properties()[3]  = particle_density;
  pit1->get_properties()[4]  = 0.01;
  pit1->get_properties()[5]  = 0;
  pit1->get_properties()[6]  = 0;
  pit1->get_properties()[7]  = 0;
  pit1->get_properties()[8]  = 0;
  pit1->get_properties()[9]  = 0;
  pit1->get_properties()[10] = 0;
  pit1->get_properties()[11] = 0;
  pit1->get_properties()[12] = 0;
  pit1->get_properties()[13] = 0;
  pit1->get_properties()[14] = 0;
  pit1->get_properties()[15] = 0;
  pit1->get_properties()[16] = 1;
  pit1->get_properties()[17] = 1;

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[0]  = id2;
  pit2->get_properties()[1]  = 1;
  pit2->get_properties()[2]  = particle_diameter;
  pit2->get_properties()[3]  = particle_density;
  pit2->get_properties()[4]  = 0;
  pit2->get_properties()[5]  = 0;
  pit2->get_properties()[6]  = 0;
  pit2->get_properties()[7]  = 0;
  pit2->get_properties()[8]  = 0;
  pit2->get_properties()[9]  = 0;
  pit2->get_properties()[10] = 0;
  pit2->get_properties()[11] = 0;
  pit2->get_properties()[12] = 0;
  pit2->get_properties()[13] = 0;
  pit2->get_properties()[14] = 0;
  pit2->get_properties()[15] = 0;
  pit2->get_properties()[16] = 1;
  pit2->get_properties()[17] = 1;

  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
    pairs;
  ppbs.find_PP_Contact_Pairs(particle_handler, cellNeighbor, pairs);

  std::vector<std::map<int, pp_contact_info_struct<dim>>> inContactInfo(
    num_Particles);

  ppfs.pp_Fine_Search(pairs, inContactInfo, dt);

  PPNonLinearForce<dim> ppnlf;
  ppnlf.calculate_pp_contact_force(&inContactInfo, dem_parameters);

  auto particle = particle_handler.begin();

  deallog << "The contact force vector for particle 1 is: "
          << particle->get_properties()[10] << " "
          << particle->get_properties()[11] << " "
          << particle->get_properties()[12] << " N " << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
