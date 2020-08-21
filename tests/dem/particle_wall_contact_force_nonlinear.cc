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

// In this test, the performance of non-linear (Hertzian) particle-wall contact
// force  is checked

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
#include <dem/pw_contact_force.h>
#include <dem/pw_fine_search.h>
#include <dem/pw_nonlinear_force.h>

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
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  const unsigned int n_properties = 21;

  // Defining general simulation parameters
  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt                                          = 0.00001;
  double         particle_diameter                           = 0.005;
  int            particle_density                            = 2500;
  dem_parameters.physical_properties.Youngs_modulus_particle = 50000000;
  dem_parameters.physical_properties.Youngs_modulus_wall     = 50000000;
  dem_parameters.physical_properties.Poisson_ratio_particle  = 0.3;
  dem_parameters.physical_properties.Poisson_ratio_wall      = 0.3;
  dem_parameters.physical_properties.restitution_coefficient_particle = 0.5;
  dem_parameters.physical_properties.restitution_coefficient_wall     = 0.5;
  dem_parameters.physical_properties.friction_coefficient_particle    = 0.5;
  dem_parameters.physical_properties.friction_coefficient_wall        = 0.5;
  dem_parameters.physical_properties.rolling_friction_particle        = 0.1;
  dem_parameters.physical_properties.rolling_friction_wall            = 0.1;

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

  // Finding boundary cells
  std::vector<typename Triangulation<dim>::active_cell_iterator>
                                                 boundary_cells_with_faces;
  std::map<int, boundary_cells_info_struct<dim>> boundary_cell_information;
  FindBoundaryCellsInformation<dim>              boundary_cells_object;
  boundary_cell_information =
    boundary_cells_object.find_boundary_cells_information(
      boundary_cells_with_faces, tr);

  // Calling broad search
  PWBroadSearch<dim> broad_search_object;
  std::unordered_map<
    int,
    std::unordered_map<
      int,
      std::tuple<Particles::ParticleIterator<dim>, Tensor<1, dim>, Point<dim>>>>
    pw_contact_list;
  broad_search_object.find_PW_Contact_Pairs(boundary_cell_information,
                                            particle_handler,
                                            pw_contact_list);

  // Calling fine search
  PWFineSearch<dim> fine_search_object;
  std::map<int, std::map<int, pw_contact_info_struct<dim>>>
    pw_contact_information;
  fine_search_object.pw_Fine_Search(pw_contact_list, pw_contact_information);

  // Calling non-linear force
  PWNonLinearForce<dim> force_object;
  force_object.calculate_pw_contact_force(&pw_contact_information,
                                          dem_parameters,
                                          dt);

  // Output
  auto particle = particle_handler.begin();
  deallog << "The contact force acting on particle 1 is: "
          << particle->get_properties()[DEM::PropertiesIndex::force_x] << " N "
          << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
