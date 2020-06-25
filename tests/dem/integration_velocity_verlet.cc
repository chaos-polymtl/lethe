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

// This test checks the performance of the velocity verlet integrator class

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <dem/velocity_verlet_integrator.h>

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

  // Defining simulation general parameters
  const unsigned int n_properties = 24;
  Tensor<1, dim>     g{{0, 0, -9.81}};
  double             dt = 0.00001;

  // Defning particle handler
  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  // inserting one particle at x = 0 , y = 0 and z = 0 m
  // initial velocity of particles = 0, 0, 0 m/s
  // gravitational acceleration = 0, 0, -9.81 m/s2
  Point<3> position1 = {0, 0, 0};
  int      id        = 0;

  Particles::Particle<dim> particle1(position1, position1, id);
  typename Triangulation<dim>::active_cell_iterator particle_cell =
    GridTools::find_active_cell_around_point(tr, particle1.get_location());

  // Inserting one particle and defining its properties
  Particles::ParticleIterator<dim> pit =
    particle_handler.insert_particle(particle1, particle_cell);
  pit->get_properties()[0] = id;
  pit->get_properties()[1] = 1;
  pit->get_properties()[2] = 0.005;
  pit->get_properties()[3] = 2500;
  // Velocity
  pit->get_properties()[4] = 0;
  pit->get_properties()[5] = 0;
  pit->get_properties()[6] = 0;
  // Acceleration
  pit->get_properties()[7] = 0;
  pit->get_properties()[8] = 0;
  pit->get_properties()[9] = -9.81;
  // Force
  pit->get_properties()[10] = 0;
  pit->get_properties()[11] = 0;
  pit->get_properties()[12] = 0;
  // w
  pit->get_properties()[13] = 0;
  pit->get_properties()[14] = 0;
  pit->get_properties()[15] = 0;
  // mass and moment of inertia
  pit->get_properties()[16] = 1;
  pit->get_properties()[17] = 1;

  // Calling velocity verlet integrator
  VelocityVerletIntegrator<dim> integration_object;
  integration_object.integrate(particle_handler, g, dt);

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
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
