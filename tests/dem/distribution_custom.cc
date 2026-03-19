// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserting particles following a custom distribution. At the end, the
 * cumulative density function is reconstructed and should fall back to the
 * input probability function used.
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
          DistributionWeightingType weighting_type,
          ProbabilityFunctionType   function_type,
          bool                      interpolate>
void
test(const std::vector<double> &diameter_list,
     const std::vector<double> &probability_list)
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

  InsertionInfo<dim> &insert_info = dem_parameters.insertion_info;

  LagrangianPhysicalProperties &lpp =
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
  lpp.particle_custom_diameter.push_back(diameter_list);
  lpp.particle_custom_probability.push_back(probability_list);
  lpp.seed_for_distributions.push_back(10);
  lpp.diameter_min_cutoff.push_back(-1.);
  lpp.diameter_max_cutoff.push_back(-1.);
  lpp.distribution_weighting_type.push_back(weighting_type);
  lpp.custom_probability_function_type.push_back(function_type);
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
    interpolate));

  // Calling volume insertion
  InsertionVolume<dim, PropertiesIndex> insertion_object(
    distribution_object_container,
    tr,
    dem_parameters,
    distribution_object_container[0]->find_max_diameter());

  insertion_object.insert(particle_handler, tr, dem_parameters);

  // Output
  if constexpr (weighting_type == DistributionWeightingType::number_based)
    deallog << "Numbered weighted custom distribution " << std::endl;
  if constexpr (weighting_type == DistributionWeightingType::volume_based)
    deallog << "Volume weighted custom distribution " << std::endl;

  unsigned int total_particle_number = particle_handler.n_global_particles();

  double           total_volume = 0.;
  constexpr double one_over_6   = 1. / 6.;

  std::vector<double> number_based_cdf(diameter_list.size());
  std::vector<double> volume_based_cdf(diameter_list.size());

  // Compute the total volume of particle
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      auto         particle_properties = particle->get_properties();
      const double dp = particle_properties[PropertiesIndex::dp];
      total_volume += M_PI * Utilities::fixed_power<3>(dp) * one_over_6;
    }

  // Loop over every particle, check if the diameter is smaller or equal to the
  // input diameter value, print the CDF.
  for (const double diameter_value_i : diameter_list)
    {
      double n_smaller_or_equal_particles                  = 0;
      double volume_occupied_by_smaller_or_equal_particles = 0.;
      for (auto particle = particle_handler.begin();
           particle != particle_handler.end();
           ++particle)
        {
          auto         particle_properties = particle->get_properties();
          const double dp = particle_properties[PropertiesIndex::dp];
          if (dp <= diameter_value_i)
            {
              n_smaller_or_equal_particles += 1.;
              volume_occupied_by_smaller_or_equal_particles +=
                M_PI * Utilities::fixed_power<3>(dp) * one_over_6;
            }
        }
      // Output
      if constexpr (weighting_type == DistributionWeightingType::number_based)
        {
          const double number_fraction =
            n_smaller_or_equal_particles / total_particle_number;

          deallog << "Number fraction smaller then " << diameter_value_i
                  << " : " << number_fraction << std::endl;
        }

      if constexpr (weighting_type == DistributionWeightingType::volume_based)
        deallog << "Volume fraction smaller then " << diameter_value_i << " : "
                << volume_occupied_by_smaller_or_equal_particles / total_volume
                << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      // Discrete
      const std::vector<double> d_list_1 = {0.0025, 0.0050};
      const std::vector<double> p_list_1 = {0.5, 0.5};
      initlog();
      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::volume_based,
           ProbabilityFunctionType::PDF,
           false>(d_list_1, p_list_1);
      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::number_based,
           ProbabilityFunctionType::PDF,
           false>(d_list_1, p_list_1);

      // Interpolate
      const std::vector<double> d_list_2 = {
        0.0025, 0.0030, 0.0040, 0.0050, 0.0055, 0.0059};
      const std::vector<double> p_list_2 = {0.1, 0.2, 0.35, 0.60, 0.76, 1.0};
      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::volume_based,
           ProbabilityFunctionType::CDF,
           false>(d_list_2, p_list_2);

      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::number_based,
           ProbabilityFunctionType::CDF,
           false>(d_list_2, p_list_2);
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
