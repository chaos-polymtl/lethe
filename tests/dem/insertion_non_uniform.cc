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
 * @brief Inserting one particle using uniform insertion class.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>

// Lethe
#include <dem/dem_solver_parameters.h>
#include <dem/non_uniform_insertion.h>

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

  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  // Defining simulation general parameters
  dem_parameters.insertion_info.x_min                             = -0.05;
  dem_parameters.insertion_info.y_min                             = -0.05;
  dem_parameters.insertion_info.z_min                             = -0.05;
  dem_parameters.insertion_info.x_max                             = 0.05;
  dem_parameters.insertion_info.y_max                             = 0.05;
  dem_parameters.insertion_info.z_max                             = 0.05;
  dem_parameters.insertion_info.inserted_this_step                = 10;
  dem_parameters.insertion_info.distance_threshold                = 2;
  dem_parameters.physical_properties.particle_type_number         = 1;
  dem_parameters.physical_properties.particle_average_diameter[0] = 0.005;
  dem_parameters.physical_properties.particle_size_std[0]         = 0;
  dem_parameters.physical_properties.density[0]                   = 2500;
  dem_parameters.physical_properties.number[0]                    = 10;
  dem_parameters.insertion_info.random_number_range               = 0.75;
  dem_parameters.insertion_info.random_number_seed                = 19;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // Calling uniform insertion
  NonUniformInsertion<dim> insertion_object(
    dem_parameters,
    dem_parameters.physical_properties.particle_average_diameter[0]);
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
