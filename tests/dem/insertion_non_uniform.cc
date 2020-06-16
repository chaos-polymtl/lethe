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

// Inserting one particle using uniform insertion class

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>

#include <dem/dem_solver_parameters.h>
#include <dem/non_uniform_insertion.h>

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

  // Defining simulation general parameters
  const unsigned int n_properties                        = 21;
  dem_parameters.insertionInfo.x_min                     = -0.05;
  dem_parameters.insertionInfo.y_min                     = -0.05;
  dem_parameters.insertionInfo.z_min                     = -0.05;
  dem_parameters.insertionInfo.x_max                     = 0.05;
  dem_parameters.insertionInfo.y_max                     = 0.05;
  dem_parameters.insertionInfo.z_max                     = 0.05;
  dem_parameters.insertionInfo.inserted_this_step        = 10;
  dem_parameters.insertionInfo.distance_threshold        = 2;
  dem_parameters.physicalProperties.diameter             = 0.005;
  dem_parameters.physicalProperties.density              = 2500;
  dem_parameters.simulationControl.total_particle_number = 10;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(tr, mapping, n_properties);

  // Calling uniform insertion
  NonUniformInsertion<dim> insertion_object(dem_parameters);
  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Output
  int particle_number = 1;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      deallog << "Particle " << particle_number
              << " is inserted at: " << particle->get_location()[0] << " "
              << particle->get_location()[1] << " "
              << particle->get_location()[2] << " " << std::endl;
    }
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();
  test<3>();
}
