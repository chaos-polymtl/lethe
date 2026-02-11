// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Insertion of particles using the plane insertion class.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

// Lethe
#include <dem/dem_solver_parameters.h>
#include <dem/insertion_plane.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 2;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  tr.refine_global(refinement_number);

  MappingQ<dim>            mapping(1);
  DEMSolverParameters<dim> dem_parameters;

  InsertionInfo<dim>           &insert_info = dem_parameters.insertion_info;
  LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;

  // Defining simulation general parameters
  // Insertion info
  insert_info.insertion_plane_normal_vector = Tensor<1, 3>({0., 1., 0.});
  insert_info.insertion_plane_point         = Point<3>({0., 1.75, 0});
  insert_info.distance_threshold            = 0.25;
  insert_info.insertion_maximum_offset      = 0.2;
  insert_info.seed_for_insertion            = 19;
  insert_info.insertion_frequency           = 2;

  // Lagrangian physical properties
  lpp.particle_type_number = 1;
  lpp.particle_average_diameter.push_back(0.2);
  lpp.distribution_type.push_back(SizeDistributionType::uniform);
  lpp.particle_size_std.push_back(0);
  lpp.density_particle.push_back(2500);
  lpp.number.push_back(16);

  // Calling uniform insertion
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;
  distribution_object_container.push_back(std::make_shared<UniformDistribution>(
    dem_parameters.lagrangian_physical_properties
      .particle_average_diameter[0]));

  // Calling plane insertion
  InsertionPlane<dim, PropertiesIndex> insertion_object(
    distribution_object_container, tr, dem_parameters);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

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
      test<3, DEM::DEMProperties::PropertiesIndex>();
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
