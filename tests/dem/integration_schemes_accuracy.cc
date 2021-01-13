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
  Tensor<1, dim>     g{{0, 0, 0}};
  double             dt1             = 0.1;
  double             dt2             = 0.05;
  const unsigned int time_step_ratio = dt1 / dt2;

  // Defning particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  int    id              = 0;
  double particle_mass   = 1;
  double spring_constant = 1;

  // Initial condition
  double   t         = 0;
  double   x0        = 0.3;
  Point<3> position1 = {0, 0, x0};
  double   t_final   = 0.999999999;
  double   x_analytical;
  double   particle_axial_position_error_Euler_dt1;
  double   particle_axial_position_error_Euler_dt2;
  double   particle_axial_position_error_Verlet_dt1;
  double   particle_axial_position_error_Verlet_dt2;
  double   particle_axial_position_error_Gear3_dt1;
  double   particle_axial_position_error_Gear3_dt2;

  // Output Analytical
  x_analytical = x0 * cos(sqrt(spring_constant / particle_mass) * t_final);

  Particles::Particle<dim> particle0(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle0_cell =
    GridTools::find_active_cell_around_point(tr, particle0.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit0 =
    particle_handler.insert_particle(particle0, particle0_cell);

  pit0->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit0->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit0->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit0->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit0->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit0->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit0->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Calling integrators
  ExplicitEulerIntegrator<dim>  explicit_euler_object;
  VelocityVerletIntegrator<dim> velocity_verlet_object;
  Gear3Integrator<dim>          gear3_integration_object;


  // Explicit Euler
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)

    {
      auto particle_properties = particle_iterator->get_properties();

      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          explicit_euler_object.integrate_pre_force(particle_handler, g, dt1);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          explicit_euler_object.integrate_post_force(particle_handler, g, dt1);

          t += dt1;
        }
      particle_axial_position_error_Euler_dt1 =
        particle_iterator->get_location()[2] - x_analytical;
    }

  particle_handler.clear_particles();
  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle1_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit1 =
    particle_handler.insert_particle(particle1, particle1_cell);


  pit1->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit1->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit1->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();

      t = 0;
      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          explicit_euler_object.integrate_pre_force(particle_handler, g, dt2);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          explicit_euler_object.integrate_post_force(particle_handler, g, dt2);

          t += dt2;
        }
      particle_axial_position_error_Euler_dt2 =
        particle_iterator->get_location()[2] - x_analytical;
    }
  deallog << "Explicit Euler is a "
          << (particle_axial_position_error_Euler_dt1 /
              particle_axial_position_error_Euler_dt2) /
               time_step_ratio
          << " order integration scheme" << std::endl;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle2(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle2_cell =
    GridTools::find_active_cell_around_point(tr, particle2.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit2 =
    particle_handler.insert_particle(particle2, particle2_cell);

  pit2->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit2->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit2->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit2->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit2->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Output Velocity Verlet
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();
      t                        = 0;
      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          velocity_verlet_object.integrate_pre_force(particle_handler, g, dt1);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          velocity_verlet_object.integrate_post_force(particle_handler, g, dt1);

          t += dt1;
        }
      particle_axial_position_error_Verlet_dt1 =
        particle_iterator->get_location()[2] - x_analytical;
    }

  particle_handler.clear_particles();
  Particles::Particle<dim> particle3(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle3_cell =
    GridTools::find_active_cell_around_point(tr, particle3.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit3 =
    particle_handler.insert_particle(particle3, particle3_cell);

  pit3->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit3->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit3->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit3->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit3->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit3->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit3->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Output Velocity Verlet
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();
      t                        = 0;
      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          velocity_verlet_object.integrate_pre_force(particle_handler, g, dt2);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          velocity_verlet_object.integrate_post_force(particle_handler, g, dt2);

          t += dt2;
        }
      particle_axial_position_error_Verlet_dt2 =
        particle_iterator->get_location()[2] - x_analytical;
    }

  deallog << "Velocity Verlet is a "
          << (particle_axial_position_error_Verlet_dt1 /
              particle_axial_position_error_Verlet_dt2) /
               time_step_ratio
          << " order integration scheme" << std::endl;
  ;

  particle_handler.clear_particles();
  Particles::Particle<dim> particle4(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle4_cell =
    GridTools::find_active_cell_around_point(tr, particle4.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit4 =
    particle_handler.insert_particle(particle4, particle4_cell);

  pit4->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit4->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit4->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit4->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit4->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit4->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit4->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit4->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit4->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit4->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit4->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Output Gear3
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();
      t                        = 0;
      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          gear3_integration_object.integrate_pre_force(particle_handler,
                                                       g,
                                                       dt1);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          gear3_integration_object.integrate_post_force(particle_handler,
                                                        g,
                                                        dt1);

          t += dt1;
        }
      particle_axial_position_error_Gear3_dt1 =
        particle_iterator->get_location()[2] - x_analytical;
    }

  particle_handler.clear_particles();
  Particles::Particle<dim> particle5(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle5_cell =
    GridTools::find_active_cell_around_point(tr, particle5.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit5 =
    particle_handler.insert_particle(particle5, particle5_cell);

  pit5->get_properties()[DEM::PropertiesIndex::v_x]         = 0;
  pit5->get_properties()[DEM::PropertiesIndex::v_y]         = 0;
  pit5->get_properties()[DEM::PropertiesIndex::v_z]         = 0;
  pit5->get_properties()[DEM::PropertiesIndex::acc_x]       = 0;
  pit5->get_properties()[DEM::PropertiesIndex::acc_y]       = 0;
  pit5->get_properties()[DEM::PropertiesIndex::acc_z]       = 0;
  pit5->get_properties()[DEM::PropertiesIndex::force_x]     = 0;
  pit5->get_properties()[DEM::PropertiesIndex::force_y]     = 0;
  pit5->get_properties()[DEM::PropertiesIndex::force_z]     = 0;
  pit5->get_properties()[DEM::PropertiesIndex::mass]        = particle_mass;
  pit5->get_properties()[DEM::PropertiesIndex::mom_inertia] = 1;

  // Output Gear3
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      auto particle_properties = particle_iterator->get_properties();
      t                        = 0;
      while (t < t_final)
        {
          particle_properties[DEM::PropertiesIndex::acc_z] =
            -spring_constant * particle_iterator->get_location()[2];
          gear3_integration_object.integrate_pre_force(particle_handler,
                                                       g,
                                                       dt2);
          particle_properties[DEM::PropertiesIndex::force_z] =
            -spring_constant * particle_iterator->get_location()[2];
          gear3_integration_object.integrate_post_force(particle_handler,
                                                        g,
                                                        dt2);

          t += dt2;
        }
      particle_axial_position_error_Gear3_dt2 =
        particle_iterator->get_location()[2] - x_analytical;
    }

  deallog << "Gear3 is a "
          << (particle_axial_position_error_Gear3_dt1 /
              particle_axial_position_error_Gear3_dt2) /
               time_step_ratio
          << " order integration scheme" << std::endl;
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
