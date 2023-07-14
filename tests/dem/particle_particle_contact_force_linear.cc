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
 * @brief In this test, the performance of linear (Hookean) particle-particle
 * contact force  is checked.
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

#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_fine_search.h>


// Tests (with common definitions)
#include <../tests/tests.h>

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
  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining general simulation parameters
  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt                                                  = 0.00001;
  double         particle_diameter                                   = 0.005;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_particle[0] =
    50000000;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[0] = 0.3;
  dem_parameters.lagrangian_physical_properties
    .restitution_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties
    .friction_coefficient_particle[0] = 0.5;
  dem_parameters.lagrangian_physical_properties
    .rolling_friction_coefficient_particle[0]                         = 0.1;
  dem_parameters.lagrangian_physical_properties.rolling_friction_wall = 0.1;
  dem_parameters.lagrangian_physical_properties.density_particle[0]   = 2500;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, DEM::get_number_properties());

  // Creating containers manager for finding cell neighbor and also broad and
  // fine particle-particle search objects
  DEMContactManager<dim> container_manager;

  // Finding cell neighbors
  FindCellNeighbors<dim> cell_neighbor_object;
  cell_neighbor_object.find_cell_neighbors(
    triangulation,
    container_manager.cells_local_neighbor_list,
    container_manager.cells_ghost_neighbor_list);

  // Inserting two particles in contact
  Point<3> position1 = {0.4, 0, 0};
  int      id1       = 0;
  Point<3> position2 = {0.40499, 0, 0};
  int      id2       = 1;

  Particles::Particle<dim> particle1(position1, position1, id1);
  typename Triangulation<dim>::active_cell_iterator cell1 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle1.get_location());
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

  Particles::Particle<dim> particle2(position2, position2, id2);
  typename Triangulation<dim>::active_cell_iterator cell2 =
    GridTools::find_active_cell_around_point(triangulation,
                                             particle2.get_location());
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, cell2);
  pit2->get_properties()[DEM::PropertiesIndex::type]    = 0;
  pit2->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
  pit2->get_properties()[DEM::PropertiesIndex::v_x]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::mass]    = 1;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (unsigned i = 0; i < MOI.size(); ++i)
    MOI[i] = 1;

  // Calling broad search
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      container_manager.particle_container[particle_iterator->get_id()] =
        particle_iterator;
    }

  container_manager.execute_particle_particle_broad_search(particle_handler);

  // Calling fine search
  container_manager.execute_particle_particle_fine_search(
    neighborhood_threshold);

  // Calling linear force
  ParticleParticleContactForce<
    dim,
    Parameters::Lagrangian::ParticleParticleContactForceModel::linear,
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
    linear_force_object(dem_parameters);
  linear_force_object.calculate_particle_particle_contact_force(
    container_manager, dt, torque, force);

  // Output
  auto particle = particle_handler.begin();
  deallog << "The contact force vector for particle 1 is: "
          << force[particle->get_id()][0] << " " << force[particle->get_id()][1]
          << " " << force[particle->get_id()][2] << " N " << std::endl;
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
