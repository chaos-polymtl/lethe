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
 * @brief In this test, post-collision velocity of a particle
 * is calculated after a particle-wall collision
 */

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>
#include <dem/particle_wall_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>

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

using namespace dealii;

template <int dim>
void
test(double coefficient_of_restitution)
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
  Tensor<1, dim> g{{0, 0, 0}};
  double         dt                           = 0.00000001;
  double         particle_diameter            = 0.002;
  unsigned int   rotating_wall_maximum_number = 6;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_particle[0] =
    800000000;
  dem_parameters.lagrangian_physical_properties.youngs_modulus_wall = 800000000;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[0] = 0.3;
  dem_parameters.lagrangian_physical_properties.poisson_ratio_wall        = 0.3;
  dem_parameters.lagrangian_physical_properties
    .restitution_coefficient_particle[0] = coefficient_of_restitution;
  dem_parameters.lagrangian_physical_properties.restitution_coefficient_wall =
    coefficient_of_restitution;
  dem_parameters.lagrangian_physical_properties
    .friction_coefficient_particle[0]                                     = 0.3;
  dem_parameters.lagrangian_physical_properties.friction_coefficient_wall = 0.3;
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

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // Inserting one particle in contact with wall
  Point<dim>               position = {-0.999, 0, 0};
  int                      id       = 0;
  Particles::Particle<dim> particle(position, position, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle.get_location());
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle, particle_cell);
  pit1->get_properties()[DEM::PropertiesIndex::type]    = 0;
  pit1->get_properties()[DEM::PropertiesIndex::dp]      = particle_diameter;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]     = -0.1;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass] =
    M_PI * particle_diameter * particle_diameter * particle_diameter / 6;


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

  // P-W broad search
  ParticleWallBroadSearch<dim> particle_wall_broad_search_object;
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  particle_wall_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // P-W fine search
  ParticleWallFineSearch<dim> particle_wall_fine_search_object;
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                                  particle_wall_contact_information;
  ParticleWallNonLinearForce<dim> particle_wall_force_object(
    dem_parameters.boundary_conditions.boundary_translational_velocity,
    dem_parameters.boundary_conditions.boundary_rotational_speed,
    dem_parameters.boundary_conditions.boundary_rotational_vector,
    grid_radius,
    dem_parameters);
  VelocityVerletIntegrator<dim> integrator_object;

  auto particle1 = particle_handler.begin();

  for (double time = 0; time < 0.0001; time += dt)
    {
      // If particle and wall are in contact
      particle_wall_fine_search_object.particle_wall_fine_search(
        particle_wall_contact_list, particle_wall_contact_information);

      particle_wall_force_object.calculate_particle_wall_contact_force(
        particle_wall_contact_information, dt, torque, force);
      integrator_object.integrate(particle_handler, g, dt, torque, force, MOI);
    }

  deallog << "Coefficient of restitution is " << coefficient_of_restitution
          << " and the velocity of particle "
             "before collision is 0.1, the velocity of particle after "
             "collision is: "
          << particle1->get_properties()[DEM::PropertiesIndex::v_x]
          << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>(0.9);
  test<3>(1);
}
