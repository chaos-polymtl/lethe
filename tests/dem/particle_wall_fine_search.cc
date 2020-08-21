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

// In this test, the performance of particle-wall fine search is investigated

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
#include <dem/find_boundary_cells_information.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_info_struct.h>
#include <dem/pw_fine_search.h>

#include <iostream>
#include <vector>

#include "../tests.h"

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  double             particle_diameter = 0.005;
  int                particle_density  = 2500;
  const unsigned int n_properties      = 21;
  double             dt                = 0.00001;

  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  // Inserting one particle in contact with a wall
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

  // Calling find_boundary_cells_information function to find the information of
  // boundary cells
  std::vector<typename Triangulation<dim>::active_cell_iterator>
                                                 boundary_cells_with_faces;
  std::map<int, boundary_cells_info_struct<dim>> boundary_cells_information;
  FindBoundaryCellsInformation<dim>              boundary_cells_object;
  boundary_cells_information =
    boundary_cells_object.find_boundary_cells_information(
      boundary_cells_with_faces, tr);

  // Calling particle-wall broad search
  PWBroadSearch<dim> broad_search_object;
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    pw_contact_list;
  broad_search_object.find_PW_Contact_Pairs(boundary_cells_information,
                                            particle_handler,
                                            pw_contact_list);

  // Calling particle-wall fine search
  PWFineSearch<dim> fine_search_object;
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    pw_contact_information;
  fine_search_object.pw_Fine_Search(pw_contact_list, pw_contact_information);

  // Output
  for (auto pw_contact_information_iterator = pw_contact_information.begin();
       pw_contact_information_iterator != pw_contact_information.end();
       ++pw_contact_information_iterator)
    {
      auto info_iterator = pw_contact_information_iterator->second.begin();
      while (info_iterator != pw_contact_information_iterator->second.end())
        {
          deallog << "Particle " << info_iterator->second.particle->get_id()
                  << " is in contact with face " << info_iterator->first
                  << std::endl;
          deallog << "The normal vector of collision is: "
                  << info_iterator->second.normal_vector[0] << " "
                  << info_iterator->second.normal_vector[1] << " "
                  << info_iterator->second.normal_vector[2] << std::endl;
          ++info_iterator;
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
