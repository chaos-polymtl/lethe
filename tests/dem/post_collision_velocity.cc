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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_force.h>
#include <dem/pw_fine_search.h>
#include <dem/pw_nonlinear_force.h>
#include <dem/velocity_verlet_integrator.h>

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
  double         dt                                             = 0.00000001;
  double         particle_diameter                              = 0.002;
  int            particle_density                               = 2000;
  unsigned int   rotating_wall_maximum_number                   = 6;
  dem_parameters.physical_properties.particle_type_number       = 1;
  dem_parameters.physical_properties.youngs_modulus_particle[0] = 800000000;
  dem_parameters.physical_properties.youngs_modulus_wall        = 800000000;
  dem_parameters.physical_properties.poisson_ratio_particle[0]  = 0.3;
  dem_parameters.physical_properties.poisson_ratio_wall         = 0.3;
  dem_parameters.physical_properties.restitution_coefficient_particle[0] = coefficient_of_restitution;
  dem_parameters.physical_properties.restitution_coefficient_wall        = coefficient_of_restitution;
  dem_parameters.physical_properties.friction_coefficient_particle[0]    = 0.3;
  dem_parameters.physical_properties.friction_coefficient_wall           = 0.3;
  dem_parameters.physical_properties.rolling_friction_coefficient_particle[0] =
    0.1;
  dem_parameters.physical_properties.rolling_friction_wall = 0.1;

  // Initializing motion of boundaries
  Tensor<1, dim> translational_and_rotational_veclocity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      translational_and_rotational_veclocity[d] = 0;
    }
  for (unsigned int counter = 0; counter < rotating_wall_maximum_number;
       ++counter)
    {
      dem_parameters.boundary_motion.boundary_rotational_speed.insert(
        {counter, 0});
      dem_parameters.boundary_motion.boundary_translational_velocity.insert(
        {counter, translational_and_rotational_veclocity});
      dem_parameters.boundary_motion.boundary_rotational_vector.insert(
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
  pit1->get_properties()[DEM::PropertiesIndex::rho]     = particle_density;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]     = -0.1;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_x]   = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_y]   = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_z]   = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass] =
    M_PI * particle_diameter * particle_diameter * particle_diameter / 6;
  pit1->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Finding boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;

  boundary_cells_object.build(tr);

  // P-W broad search
  PWBroadSearch<dim> pw_broad_search_object;
  std::unordered_map<
    int,
    std::unordered_map<int,
                       std::tuple<Particles::ParticleIterator<dim>,
                                  Tensor<1, dim>,
                                  Point<dim>,
                                  unsigned int>>>
    pw_contact_list;
  pw_broad_search_object.find_particle_wall_contact_pairs(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    pw_contact_list);

  // P-W fine search
  PWFineSearch<dim> pw_fine_search_object;
  std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
                        pw_contact_information;
  PWNonLinearForce<dim> pw_force_object(
    dem_parameters.boundary_motion.boundary_translational_velocity,
    dem_parameters.boundary_motion.boundary_rotational_speed,
    dem_parameters.boundary_motion.boundary_rotational_vector,
    grid_radius,
    dem_parameters);
  VelocityVerletIntegrator<dim> integrator_object;

  auto particle1 = particle_handler.begin();

  for (double time = 0; time < 0.0001; time += dt)
    {
      // If particle and wall are in contact
      pw_fine_search_object.particle_wall_fine_search(pw_contact_list,
                                                      pw_contact_information);

      integrator_object.integrate_pre_force(particle_handler, g, dt);
      pw_force_object.calculate_pw_contact_force(pw_contact_information, dt);
      integrator_object.integrate_post_force(particle_handler, g, dt);
    }

  deallog << "Coefficient of restitution is " << coefficient_of_restitution << " and the velocity of particle "
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
