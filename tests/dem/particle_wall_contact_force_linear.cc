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

/**
 * @brief In this test, the performance of linear (Hookean) particle-wall
 * contact force  is checked,
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>
#include <dem/particle_wall_linear_force.h>

// Tests (with common definitions)
#include <../tests/tests.h>

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
  const double grid_radius       = 0.5 * GridTools::diameter(tr);
  int          refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt                                                  = 0.00001;
  double         particle_diameter                                   = 0.005;
  unsigned int   rotating_wall_maximum_number                        = 6;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_particle[0] =
    50000000;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_wall = 50000000;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[0] = 0.3;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_wall        = 0.3;
  dem_parameters.lagrangian_physical_properties
    .restitution_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties.restitution_coefficient_wall =
    0.5;
  dem_parameters.lagrangian_physical_properties
    .friction_coefficient_particle[0]                                     = 0.5;
  dem_parameters.lagrangian_physical_properties.friction_coefficient_wall = 0.5;
  dem_parameters.lagrangian_physical_properties
    .rolling_friction_coefficient_particle[0]                         = 0.1;
  dem_parameters.lagrangian_physical_properties.rolling_friction_wall = 0.1;
  dem_parameters.lagrangian_physical_properties.density_particle[0]   = 2500;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  // Initializing motion of boundaries
  Tensor<1, dim> translational_and_rotational_veclocity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      translational_and_rotational_veclocity[d] = 0;
    }
  for (unsigned int counter = 0; counter < rotating_wall_maximum_number;
       ++counter)
    {
      dem_parameters.boundary_conditions.boundary_rotational_speed.insert(
        {counter, 0});
      dem_parameters.boundary_conditions.boundary_translational_velocity.insert(
        {counter, translational_and_rotational_veclocity});
      dem_parameters.boundary_conditions.boundary_rotational_vector.insert(
        {counter, translational_and_rotational_veclocity});
    }

  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // Inserting one particle in contact with a wall
  Point<dim>               position1 = {-0.998, 0, 0};
  int                      id1       = 0;
  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, cell1);
  pit1->get_properties()[DEM::PropertiesIndex::type]    = 0;
  pit1->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]     = 0.01;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass]    = 1;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (unsigned i = 0; i < MOI.size(); ++i)
    MOI[i] = 1;

  // Finding boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::vector<unsigned int>     outlet_boundaries;
  boundary_cells_object.build(
    tr,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Calling broad search
  ParticleWallBroadSearch<dim> broad_search_object;
  // P-W broad search
  ParticleWallBroadSearch<dim> particle_wall_broad_search_object;
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);
  broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // Calling fine search
  ParticleWallFineSearch<dim> fine_search_object;
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_contact_information;
  fine_search_object.particle_wall_fine_search(
    particle_wall_contact_list, particle_wall_contact_information);

  // Calling linear force
  ParticleWallLinearForce<dim> force_object(
    dem_parameters.boundary_conditions.boundary_translational_velocity,
    dem_parameters.boundary_conditions.boundary_rotational_speed,
    dem_parameters.boundary_conditions.boundary_rotational_vector,
    grid_radius,
    dem_parameters);
  force_object.calculate_particle_wall_contact_force(
    particle_wall_contact_information, dt, torque, force);

  // Output
  auto particle = particle_handler.begin();
  deallog << "The contact force acting on particle 1 is: "
          << force[particle->get_id()][0] << " N " << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<3>();
    }
  catch (std::exception &exc)
    {
      try
        {
          std::cerr << std::endl
                    << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "Exception on processing: " << std::endl
                    << exc.what() << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          return 1;
        }
      catch (...)
        {
          std::cerr << std::endl
                    << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          std::cerr << "Unknown exception!" << std::endl
                    << "Aborting!" << std::endl
                    << "----------------------------------------------------"
                    << std::endl;
          return 1;
        }
      return 0;
    }
}
