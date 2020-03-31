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
#include <dem/pw_nonlinear_force.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/find_boundary_cells_information.h"
#include "dem/pw_broad_search.h"
#include "dem/pw_contact_force.h"
#include "dem/pw_fine_search.h"
#include "dem/velocity_verlet_integrator.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  int                num_particles = 1;
  const unsigned int n_properties  = 21;

  Tensor<1, dim> g{{0, 0, 0}};
  double         dt = 0.000001;

  double particle_diameter                                  = 0.001;
  int    particle_density                                   = 7850;
  dem_parameters.physicalProperties.Youngs_modulus_particle = 200000000000;
  dem_parameters.physicalProperties.Youngs_modulus_wall     = 200000000000;
  dem_parameters.physicalProperties.Poisson_ratio_particle  = 0.3;
  dem_parameters.physicalProperties.Poisson_ratio_wall      = 0.3;
  dem_parameters.physicalProperties.restitution_coefficient_particle = 0.5;
  dem_parameters.physicalProperties.restitution_coefficient_wall     = 0.5;
  dem_parameters.physicalProperties.friction_coefficient_particle    = 0.3;
  dem_parameters.physicalProperties.friction_coefficient_wall        = 0.3;
  dem_parameters.physicalProperties.rolling_friction_particle        = 0.1;
  dem_parameters.physicalProperties.rolling_friction_wall            = 0.1;

  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  Point<dim> position1 = {-0.999, 0, 0};
  int        id1       = 0;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[0]  = id1;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = particle_diameter;
  pit1->get_properties()[3]  = particle_density;
  pit1->get_properties()[4]  = -1.0;
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

  std::vector<boundary_cells_info_struct<dim>> boundaryCellInfo;
  FindBoundaryCellsInformation<dim>            boundary_cells_object;
  boundaryCellInfo = boundary_cells_object.find_boundary_cells_information(tr);

  PWBroadSearch<dim> pw1;
  std::vector<std::tuple<std::pair<Particles::ParticleIterator<dim>, int>,
                         Tensor<1, dim>,
                         Point<dim>>>
    pwContactList(num_particles);

  pw1.find_PW_Contact_Pairs(boundaryCellInfo, particle_handler, pwContactList);

  PWFineSearch<dim>                                       pw2;
  std::vector<std::map<int, pw_contact_info_struct<dim>>> pwContactInfo(
    num_particles);
  PWNonLinearForce<dim>         pwcf1;
  VelocityVerletIntegrator<dim> Integ1;

  for (double time = 0; time < 0.00115; time += dt)
    {
      auto particle = particle_handler.begin();
      particle->get_properties()[DEM::PropertiesIndex::force_x] = 0;
      particle->get_properties()[DEM::PropertiesIndex::force_y] = 0;
      particle->get_properties()[DEM::PropertiesIndex::force_z] = 0;

      pw2.pw_Fine_Search(pwContactList, pwContactInfo, dt);
      auto pw_pairs_in_contact_iterator = pwContactInfo.begin();
      auto pw_contact_information_iterator =
        pw_pairs_in_contact_iterator->begin();

      pwcf1.calculate_pw_contact_force(&pwContactInfo, dem_parameters);

      Integ1.integrate(particle_handler, g, dt);
      if (pw_contact_information_iterator->second.normal_overlap > 0)
        {
          deallog << " "
                  << pw_contact_information_iterator->second.normal_overlap
                  << " "
                  << particle->get_properties()[DEM::PropertiesIndex::force_x]
                  << std::endl;
        }
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
