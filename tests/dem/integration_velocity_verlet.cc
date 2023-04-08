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
 * @brief This test checks the performance of the velocity verlet integrator
 * class.
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
#include <core/dem_properties.h>

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
  Tensor<1, dim> g{{0, 0, -9.81}};
  double         dt = 0.00001;

  // Defning particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // inserting one particle at x = 0 , y = 0 and z = 0 m
  // initial velocity of particles = 0, 0, 0 m/s
  // gravitational acceleration = 0, 0, -9.81 m/s2
  Point<3> position1 = {0, 0, 0};
  int      id        = 0;

  DEMSolverParameters<dim> dem_parameters;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.density_particle[0]  = 2500;

  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit =
    particle_handler.insert_particle(particle1, particle_cell);

  pit->get_properties()[DEM::PropertiesIndex::type]    = 1;
  pit->get_properties()[DEM::PropertiesIndex::dp]      = 0.005;
  pit->get_properties()[DEM::PropertiesIndex::v_x]     = 0.;
  pit->get_properties()[DEM::PropertiesIndex::v_y]     = 0.;
  pit->get_properties()[DEM::PropertiesIndex::v_z]     = 0.;
  pit->get_properties()[DEM::PropertiesIndex::omega_x] = 0.;
  pit->get_properties()[DEM::PropertiesIndex::omega_y] = 0.;
  pit->get_properties()[DEM::PropertiesIndex::omega_z] = 0.;
  pit->get_properties()[DEM::PropertiesIndex::mass]    = 1.;

  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;
  torque.push_back(Tensor<1, dim>({0, 0, 0}));
  force.push_back(Tensor<1, dim>({0, 0, 0}));
  MOI.push_back(1);

  // Calling velocity verlet integrator
  VelocityVerletIntegrator<dim> integration_object;
  integration_object.integrate(particle_handler, g, dt, torque, force, MOI);

  // Output
  for (auto particle_iterator = particle_handler.begin();
       particle_iterator != particle_handler.end();
       ++particle_iterator)
    {
      deallog << "The new position of the particle in z direction after " << dt
              << " seconds is: " << particle_iterator->get_location()[2]
              << std::endl;
    }
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
