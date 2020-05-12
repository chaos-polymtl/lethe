/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

// In this test, the performance of particle-particle fine search class is
// evaluated

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
#include <dem/find_cell_neighbors.h>
#include <dem/pp_broad_search.h>
#include <dem/pp_contact_info_struct.h>
#include <dem/pp_fine_search.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  const unsigned int              n_properties = 21;
  Particles::ParticleHandler<dim> particle_handler(triangulation,
                                                   mapping,
                                                   n_properties);

  // Finding cell neighbors
  std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
    cell_neighbors_list;

  FindCellNeighbors<dim> cn1;
  cell_neighbors_list = cn1.find_cell_neighbors(triangulation);

  // Creating broad and fine particle-particle search objects
  PPBroadSearch<dim> broad_search_object;
  PPFineSearch<dim>  fine_search_obejct;

  // Inserting two particles in contact
  Point<3> position1 = {0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.40499, 0, 0};
  int      id2       = 1;
  float    dt        = 0.00001;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[0]  = id1;
  pit1->get_properties()[1]  = 1;
  pit1->get_properties()[2]  = 0.005;
  pit1->get_properties()[3]  = 2500;
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

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[0]  = id2;
  pit2->get_properties()[1]  = 1;
  pit2->get_properties()[2]  = 0.005;
  pit2->get_properties()[3]  = 2500;
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

  // Calling broad search
  std::vector<std::pair<Particles::ParticleIterator<dim>,
                        Particles::ParticleIterator<dim>>>
    pairs;
  broad_search_object.find_PP_Contact_Pairs(particle_handler,
                                            cell_neighbors_list,
                                            pairs);

  // Calling fine search
  std::map<int, std::map<int, pp_contact_info_struct<dim>>>
    contact_pairs_information;

  fine_search_obejct.pp_Fine_Search(pairs, contact_pairs_information, dt);

  // Output
  for (auto in_contact_info_iterator = contact_pairs_information.begin();
       in_contact_info_iterator != contact_pairs_information.end();
       ++in_contact_info_iterator)
    {
      auto information_iterator = (&in_contact_info_iterator->second)->begin();
      while (information_iterator != (&in_contact_info_iterator->second)->end())
        {
          deallog << "The normal overlap of contacting paritlce is: "
                  << information_iterator->second.normal_overlap << std::endl;
          deallog << "The normal vector of collision is: "
                  << information_iterator->second.normal_vector[0] << " "
                  << information_iterator->second.normal_vector[1] << " "
                  << information_iterator->second.normal_vector[2] << std::endl;
          deallog << "Normal relative velocity at contact point is: "
                  << information_iterator->second.normal_relative_velocity
                  << std::endl;
          deallog << "The tangential overlap of contacting paritlce is: "
                  << information_iterator->second.tangential_overlap[0] << " "
                  << information_iterator->second.tangential_overlap[1] << " "
                  << information_iterator->second.tangential_overlap[2]
                  << std::endl;
          deallog
            << "Tangential relative velocity at contact point is: "
            << information_iterator->second.tangential_relative_velocity[0]
            << " "
            << information_iterator->second.tangential_relative_velocity[1]
            << " "
            << information_iterator->second.tangential_relative_velocity[2]
            << std::endl;
          ++information_iterator;
        }
    }
}

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  test<3>();
}
