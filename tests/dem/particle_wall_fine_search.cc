#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

#include "../tests.h"
#include "dem/dem_solver_parameters.h"
#include "dem/find_boundary_cells_information.h"
#include "dem/pw_broad_search.h"
#include "dem/pw_contact_info_struct.h"
#include "dem/pw_fine_search.h"

using namespace dealii;

template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tr, -1, 1, true);
  int numRef = 2;
  tr.refine_global(numRef);
  int           cellNum = tr.n_active_cells();
  MappingQ<dim> mapping(1);

  double             particle_diameter = 0.005;
  int                particle_density  = 2500;
  const unsigned int n_properties      = 21;
  double             dt                = 0.00001;

  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  Point<dim>               position1 = {-0.998, 0, 0};
  int                      id1       = 0;
  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[0]  = id1;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = particle_diameter;
  pit1->get_properties()[3]  = particle_density;
  pit1->get_properties()[4]  = 0;
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
  std::map<int,
           std::tuple<std::pair<typename Particles::ParticleIterator<dim>, int>,
                      Tensor<1, dim>,
                      Point<dim>>>
    pwContactList;

  pw1.find_PW_Contact_Pairs(boundaryCellInfo, particle_handler, pwContactList);

  PWFineSearch<dim> pw2;

  std::map<int, std::map<int, pw_contact_info_struct<dim>>> pwContactInfo;

  pw2.pw_Fine_Search(pwContactList, pwContactInfo, dt);
  for (auto pwContactInfo_iterator = pwContactInfo.begin();
       pwContactInfo_iterator != pwContactInfo.end();
       ++pwContactInfo_iterator)
    {
      auto info_it = pwContactInfo_iterator->second.begin();
      while (info_it != pwContactInfo_iterator->second.end())
        {
          deallog << "The normal overlap of contacting paritlce-wall is: "
                  << info_it->second.normal_overlap << std::endl;
          deallog << "The normal vector of collision is: "
                  << info_it->second.normal_vector[0] << " "
                  << info_it->second.normal_vector[1] << " "
                  << info_it->second.normal_vector[2] << std::endl;
          deallog << "Normal relative velocity at contact point is: "
                  << info_it->second.normal_relative_velocity << std::endl;
          deallog << "The tangential overlap of contacting paritlce-wall is: "
                  << info_it->second.tangential_overlap[0] << " "
                  << info_it->second.tangential_overlap[1] << " "
                  << info_it->second.tangential_overlap[2] << std::endl;
          deallog << "Tangential relative velocity at contact point is: "
                  << info_it->second.tangential_relative_velocity[0] << " "
                  << info_it->second.tangential_relative_velocity[1] << " "
                  << info_it->second.tangential_relative_velocity[2]
                  << std::endl;
          ++info_it;
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
