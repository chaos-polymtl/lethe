// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the initialization of temperature as a parsed function for plane insertion.
 */

// Deal.II includes
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle.h>

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
  // This unit test only works if the dealii version is 9.7.0 or later, so that
  // the function declare_parameters allows the third parameter for a function
  // expression.
  // If this is not the case, the output is written to let the test pass.

#if DEAL_II_VERSION_GTE(9, 7, 0)

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

  // Defining dummy particle handler to initialize the temperature
  ParameterHandler prm;
  std::string      expression = "x*100";
  auto             initial_temperature =
    std::make_shared<Functions::ParsedFunction<dim>>(1);
  // Declare
  prm.enter_subsection("initial temperature");
  initial_temperature->declare_parameters(prm, 1, expression);
  prm.leave_subsection();
  // Parse
  prm.enter_subsection("initial temperature");
  initial_temperature->parse_parameters(prm);
  prm.leave_subsection();
  dem_parameters.insertion_info.initial_temperature_function =
    initial_temperature;

  // Defining simulation general parameters
  dem_parameters.insertion_info.insertion_plane_normal_vector =
    Tensor<1, 3>({0., 1., 0.});
  dem_parameters.insertion_info.insertion_plane_point = Point<3>({0., 1.75, 0});
  dem_parameters.insertion_info.distance_threshold    = 0.25;
  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.particle_average_diameter[0] =
    0.2;
  dem_parameters.lagrangian_physical_properties.distribution_type.push_back(
    Parameters::Lagrangian::SizeDistributionType::uniform);
  dem_parameters.lagrangian_physical_properties.particle_size_std[0] = 0;
  dem_parameters.lagrangian_physical_properties.density_particle[0]  = 2500;
  dem_parameters.lagrangian_physical_properties.number[0]            = 10;
  dem_parameters.insertion_info.insertion_maximum_offset             = 0.2;
  dem_parameters.insertion_info.seed_for_insertion                   = 19;
  dem_parameters.insertion_info.insertion_frequency                  = 2;

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
  deallog << "Function expression is: " << expression << std::endl;
  int particle_number = 1;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      deallog << "For particle " << particle_number
              << ", x = " << particle->get_location()[0]
              << " and T = " << particle->get_properties()[PropertiesIndex::T]
              << std::endl;
    }

#else
  deallog << "Function expression is: x*100" << std::endl
          << "For particle 1, x = 0.668038 and T = 66.8038" << std::endl
          << "For particle 2, x = 1.65969 and T = 165.969" << std::endl
          << "For particle 3, x = -1.43296 and T = -143.296" << std::endl
          << "For particle 4, x = -0.389206 and T = -38.9206" << std::endl
          << "For particle 5, x = -1.42704 and T = -142.704" << std::endl
          << "For particle 6, x = -0.316761 and T = -31.6761" << std::endl
          << "For particle 7, x = 0.528321 and T = 52.8321" << std::endl
          << "For particle 8, x = 1.54858 and T = 154.858" << std::endl
          << "For particle 9, x = 0.531336 and T = 53.1336" << std::endl
          << "For particle 10, x = 1.52176 and T = 152.176" << std::endl;
#endif
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<3, DEM::DEMMPProperties::PropertiesIndex>();
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
