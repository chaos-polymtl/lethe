// SPDX-FileCopyrightText: Copyright (c) 2020-2021, 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserting particles wit different temperatures using list insertion class.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle.h>

// Lethe
#include <../tests/dem/test_particles_functions.h>
#include <dem/dem_solver_parameters.h>
#include <dem/insertion_list.h>


template <int dim, typename PropertiesIndex>
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
  set_default_dem_parameters(1, dem_parameters);
  dem_parameters.insertion_info.list_x = {0.1, 0.2, 0.3, 0.4, 0.5};
  dem_parameters.insertion_info.list_y.resize(5);
  dem_parameters.insertion_info.list_z.resize(5);
  dem_parameters.insertion_info.list_d = {0.001, 0.002, 0.003, 0.004, 0.005};
  dem_parameters.insertion_info.list_T = {100, 200, 300, 400, 500};
  dem_parameters.insertion_info.list_vx.resize(5);
  dem_parameters.insertion_info.list_vy.resize(5);
  dem_parameters.insertion_info.list_vz.resize(5);
  dem_parameters.insertion_info.list_wx.resize(5);
  dem_parameters.insertion_info.list_wy.resize(5);
  dem_parameters.insertion_info.list_wz.resize(5);

  dem_parameters.lagrangian_physical_properties.particle_type_number = 1;
  dem_parameters.lagrangian_physical_properties.distribution_type.push_back(
    Parameters::Lagrangian::SizeDistributionType::uniform);
  dem_parameters.lagrangian_physical_properties.particle_average_diameter[0] =
    0.005;
  dem_parameters.lagrangian_physical_properties.number[0] = 5;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);
  // Calling uniform insertion
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;
  distribution_object_container.push_back(std::make_shared<UniformDistribution>(
    dem_parameters.lagrangian_physical_properties
      .particle_average_diameter[0]));

  // Calling list insertion
  InsertionList<dim, PropertiesIndex> insertion_object(
    distribution_object_container, tr, dem_parameters);

  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Output
  int particle_number = 1;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      deallog << "For particle " << particle_number
              << ", temperature is set at: "
              << particle->get_properties()[PropertiesIndex::T]
              << " K, diameter is set at: "
              << particle->get_properties()[PropertiesIndex::dp]
              << " and x is set at: " << particle->get_location()[0] << "."
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
