// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserting particles following a custom distribution.
 */

// Deal.II includes
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

// Lethe
#include <dem/dem_solver_parameters.h>
#include <dem/insertion_volume.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim,
          typename PropertiesIndex,
          DistributionWeightingType weighting_type>
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

  Parameters::Lagrangian::InsertionInfo<dim> &insert_info =
    dem_parameters.insertion_info;

  Parameters::Lagrangian::LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;

  // Defining simulation general parameters
  // Insertion info
  insert_info.insertion_box_point_1    = {-0.5, -0.5, -0.05};
  insert_info.insertion_box_point_2    = {0.5, 0.5, 0.05};
  insert_info.direction_sequence       = {0, 1, 2};
  insert_info.inserted_this_step       = 10000;
  insert_info.distance_threshold       = 2;
  insert_info.insertion_maximum_offset = 0.75;
  insert_info.seed_for_insertion       = 19;

  // Lagrangian physical properties
  lpp.particle_type_number = 1;
  lpp.distribution_type.push_back(SizeDistributionType::custom);
  lpp.particle_custom_diameter.push_back({0.0025, 0.0050});
  lpp.particle_custom_probability.push_back({0.5, 0.5});
  lpp.seed_for_distributions.push_back(10);
  lpp.diameter_min_cutoff.push_back(-1.);
  lpp.diameter_max_cutoff.push_back(-1.);
  lpp.distribution_weighting_type.push_back(
    DistributionWeightingType::volume_based);
  lpp.custom_probability_function_type.push_back(ProbabilityFunctionType::PDF);

  // lpp.diameter_min_cutoff.push_back(-1);
  // lpp.diameter_max_cutoff.push_back(-1.);
  lpp.density_particle.push_back(2500);
  lpp.number.push_back(100000);

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Calling custom distribution
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;
  distribution_object_container.push_back(std::make_shared<CustomDistribution>(
    lpp.particle_custom_diameter.at(0),
    lpp.particle_custom_probability.at(0),
    lpp.seed_for_distributions.at(0),
    lpp.diameter_min_cutoff.at(0),
    lpp.diameter_max_cutoff.at(0),
    lpp.distribution_weighting_type.at(0),
    lpp.custom_probability_function_type.at(0),
    false));

  // Calling volume insertion
  InsertionVolume<dim, PropertiesIndex> insertion_object(
    distribution_object_container,
    tr,
    dem_parameters,
    distribution_object_container[0]->find_max_diameter());

  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Output
  if constexpr (weighting_type == DistributionWeightingType::number_based)
    deallog << "Numbered weighted normal distribution " << std::endl;
  if constexpr (weighting_type == DistributionWeightingType::volume_based)
    deallog << "Volume weighted normal distribution " << std::endl;

  unsigned int     particle_number = 0;
  unsigned int     n_particle_1 = 0, n_particle_2 = 0;
  constexpr double d1 = 0.0025, d2 = 0.0050;

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      auto particle_properties = particle->get_properties();

      double dp = particle_properties[PropertiesIndex::dp];

      // deallog << "Particle " << particle_number << " diameter is: " << dp
      //         << std::endl;

      if (dp == d1)
        n_particle_1++;
      else //(dp == d2)
        n_particle_2++;
    }
  const double volume_particle_1 = n_particle_1 * std::pow(d1, 3);
  const double volume_particle_2 = n_particle_2 * std::pow(d2, 3);
  const double total_volume      = volume_particle_1 + volume_particle_2;

  deallog << "Volume fraction particle 1: " << volume_particle_1 / total_volume
          << std::endl;
  deallog << "Volume fraction particle 2: " << volume_particle_2 / total_volume
          << std::endl;
  deallog << "Number of particle 1: " << n_particle_1 << std::endl;
  deallog << "Number of particle 2: " << n_particle_2 << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      initlog();
      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::volume_based>();
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
