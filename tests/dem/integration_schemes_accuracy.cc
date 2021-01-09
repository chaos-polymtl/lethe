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
 * @brief In this test, accuracies of different integration
 * schemes are compared using a pendulum oscillation model
 */

// Deal.II includes
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

// Lethe
#include <dem/dem_properties.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/gear3_integrator.h>
#include <dem/velocity_verlet_integrator.h>

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
  int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining simulation general parameters
  Tensor<1, dim> g{{0, 0, 9.81}};
  Tensor<1, dim> g2{{0, 0, 0}};
  double         dt = 0.001;

  // Defning particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // inserting one particle1 at x = 0 , y = 0 and z = 0 m
  // initial velocity of particles = 0, 0, 0 m/s
  // gravitational acceleration = 0, 0, -9.81 m/s2
  Point<3> position1 = {0, 0, 0.3};
  int      id        = 0;

  // Initial condition
  double       teta0  = 0.3;
  const double length = 0.1;
  double       teta;

  Particles::Particle<dim> particle0(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle0_cell =
    GridTools::find_active_cell_around_point(tr, particle0.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit0 =
    particle_handler.insert_particle(particle0, particle0_cell);

  pit0->get_properties()[DEM::PropertiesIndex::type]             = 1;
  pit0->get_properties()[DEM::PropertiesIndex::dp]               = 0.005;
  pit0->get_properties()[DEM::PropertiesIndex::rho]              = 2500;
  pit0->get_properties()[DEM::PropertiesIndex::v_x]              = 0;
  pit0->get_properties()[DEM::PropertiesIndex::v_y]              = 0;
  pit0->get_properties()[DEM::PropertiesIndex::v_z]              = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_x]            = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_y]            = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_z]            = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_x]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_y]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_z]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::omega_x]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::omega_y]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::omega_z]          = 0;
  pit0->get_properties()[DEM::PropertiesIndex::mass]             = 0.001;
  pit0->get_properties()[DEM::PropertiesIndex::mom_inertia]      = 0.001;
  pit0->get_properties()[DEM::PropertiesIndex::acc_derivative_x] = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_derivative_y] = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_derivative_z] = 0;

  // Calling integrators
  ExplicitEulerIntegrator<dim>  explicit_euler_object;
  VelocityVerletIntegrator<dim> velocity_verlet_object;
  Gear3Integrator<dim>          gear3_integration_object;

  deallog << "Calculated position of particle   " << std::endl;

  // Output Analytical
  deallog << "Analytical solution:   " << std::endl;
  for (double t = 0; t < 1;)
    {
      teta = teta0 * cos(sqrt(g[2] / length) * t);

      deallog << teta << " ,";
      t += dt;
    }
  deallog << std::endl;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle1_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle1_cell);

  pit1->get_properties()[DEM::PropertiesIndex::type]             = 1;
  pit1->get_properties()[DEM::PropertiesIndex::dp]               = 0.005;
  pit1->get_properties()[DEM::PropertiesIndex::rho]              = 2500;
  pit1->get_properties()[DEM::PropertiesIndex::v_x]              = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]              = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]              = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_x]            = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_y]            = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_z]            = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_x]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_y]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_z]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_x]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_y]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::omega_z]          = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass]             = 0.001;
  pit1->get_properties()[DEM::PropertiesIndex::mom_inertia]      = 0.001;
  pit1->get_properties()[DEM::PropertiesIndex::acc_derivative_x] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_derivative_y] = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_derivative_z] = 0;

  // Output Explicit Euler
  deallog << "Explicit Euler:   " << std::endl;
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();

      for (double t = 0; t < 1;)
        {
          particle_properties[DEM::PropertiesIndex::force_z] =
            -teta0 *
            (particle_properties[DEM::PropertiesIndex::mass] * g[2] / length) *
            cos(sqrt(g[2] / length) * t);
          explicit_euler_object.integrate(particle_handler, g2, dt);

          deallog << particle_iterator->get_location()[2] << " ,";
          t += dt;
        }
    }
  deallog << std::endl;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle2(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle2_cell =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, particle2_cell);

  pit2->get_properties()[DEM::PropertiesIndex::type]             = 1;
  pit2->get_properties()[DEM::PropertiesIndex::dp]               = 0.005;
  pit2->get_properties()[DEM::PropertiesIndex::rho]              = 2500;
  pit2->get_properties()[DEM::PropertiesIndex::v_x]              = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_y]              = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_z]              = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_x]            = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_y]            = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_z]            = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_x]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_y]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_z]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_x]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_y]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::omega_z]          = 0;
  pit2->get_properties()[DEM::PropertiesIndex::mass]             = 0.001;
  pit2->get_properties()[DEM::PropertiesIndex::mom_inertia]      = 0.001;
  pit2->get_properties()[DEM::PropertiesIndex::acc_derivative_x] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_derivative_y] = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_derivative_z] = 0;

  // Output Velocity Verlet
  deallog << "Velocity Verlet:   " << std::endl;
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();

      for (double t = 0; t < 1;)
        {
          particle_properties[DEM::PropertiesIndex::force_z] =
            -teta0 *
            (particle_properties[DEM::PropertiesIndex::mass] * g[2] / length) *
            cos(sqrt(g[2] / length) * t);
          velocity_verlet_object.integrate(particle_handler, g2, dt);

          deallog << particle_iterator->get_location()[2] << " ,";
          t += dt;
        }
    }
  deallog << std::endl;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle3(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle3_cell =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, particle3_cell);

  pit3->get_properties()[DEM::PropertiesIndex::type]             = 1;
  pit3->get_properties()[DEM::PropertiesIndex::dp]               = 0.005;
  pit3->get_properties()[DEM::PropertiesIndex::rho]              = 2500;
  pit3->get_properties()[DEM::PropertiesIndex::v_x]              = 0;
  pit3->get_properties()[DEM::PropertiesIndex::v_y]              = 0;
  pit3->get_properties()[DEM::PropertiesIndex::v_z]              = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_x]            = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_y]            = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_z]            = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_x]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_y]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_z]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::omega_x]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::omega_y]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::omega_z]          = 0;
  pit3->get_properties()[DEM::PropertiesIndex::mass]             = 0.001;
  pit3->get_properties()[DEM::PropertiesIndex::mom_inertia]      = 0.001;
  pit3->get_properties()[DEM::PropertiesIndex::acc_derivative_x] = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_derivative_y] = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_derivative_z] = 0;

  // Output Gear3
  deallog << "Gear3:   " << std::endl;
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();

      for (double t = 0; t < 1;)
        {
          particle_properties[DEM::PropertiesIndex::force_z] =
            -teta0 *
            (particle_properties[DEM::PropertiesIndex::mass] * g[2] / length) *
            cos(sqrt(g[2] / length) * t);
          gear3_integration_object.integrate(particle_handler, g2, dt);

          deallog << particle_iterator->get_location()[2] << " ,";
          t += dt;
        }
    }
  deallog << std::endl;
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
