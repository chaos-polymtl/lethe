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
 */

/**
 * @brief Inserting particles following a normal distribution.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle.h>

// Lethe
#include <dem/dem_solver_parameters.h>
#include <dem/volume_insertion.h>

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
  dem_parameters.insertion_info.x_min                                = -0.5;
  dem_parameters.insertion_info.y_min                                = -0.5;
  dem_parameters.insertion_info.z_min                                = -0.05;
  dem_parameters.insertion_info.x_max                                = 0.5;
  dem_parameters.insertion_info.y_max                                = 0.5;
  dem_parameters.insertion_info.z_max                                = 0.05;
  dem_parameters.insertion_info.axis_0                               = 0;
  dem_parameters.insertion_info.axis_1                               = 1;
  dem_parameters.insertion_info.axis_2                               = 2;
  dem_parameters.insertion_info.inserted_this_step                   = 100;
  dem_parameters.insertion_info.distance_threshold                   = 2;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.distribution_type.push_back(
    Parameters::Lagrangian::SizeDistributionType::normal);
  dem_parameters.lagrangian_physical_properties.particle_average_diameter[0] =
    0.005;
  dem_parameters.lagrangian_physical_properties.particle_size_std[0] = 0.0005;
  dem_parameters.lagrangian_physical_properties.seed_for_distributions
    .push_back(10);
  dem_parameters.lagrangian_physical_properties.density_particle[0] = 2500;
  dem_parameters.lagrangian_physical_properties.number[0]           = 1000;
  dem_parameters.insertion_info.insertion_maximum_offset            = 0.75;
  dem_parameters.insertion_info.seed_for_insertion                  = 19;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, DEM::get_number_properties());

  // Calling uniform insertion
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;
  distribution_object_container.push_back(std::make_shared<NormalDistribution>(
    dem_parameters.lagrangian_physical_properties.particle_average_diameter[0],
    dem_parameters.lagrangian_physical_properties.particle_size_std[0],
    dem_parameters.lagrangian_physical_properties.seed_for_distributions[0]));

  // Calling volume insertion
  VolumeInsertion<dim> insertion_object(
    dem_parameters,
    distribution_object_container[0]->find_max_diameter(),
    distribution_object_container);

  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Output
  int particle_number = 0;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      auto particle_properties = particle->get_properties();

      double dp = particle_properties[DEM::PropertiesIndex::dp];

      deallog << "Particle " << particle_number << " diameter is: " << dp
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
