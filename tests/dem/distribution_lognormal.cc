// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Inserting particles following a normal distribution. At the end, the
 * mean and standard deviation of the inserted particles is computed. By
 * increasing the number of inserted particles, those two values should converge
 * the parameter used as inputs.
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
test(const double mu, const double sigma)
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 0;
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
  insert_info.inserted_this_step       = 1000;
  insert_info.distance_threshold       = 2;
  insert_info.insertion_maximum_offset = 0.75;
  insert_info.seed_for_insertion       = 19;

  // Lagrangian physical properties
  lpp.particle_type_number = 1;
  lpp.distribution_weighting_type.push_back(weighting_type);
  lpp.distribution_type.push_back(
    Parameters::Lagrangian::SizeDistributionType::lognormal);
  lpp.particle_average_diameter[0] = mu;
  lpp.particle_size_std[0]         = sigma;
  lpp.seed_for_distributions.push_back(10);
  lpp.diameter_min_cutoff.push_back(-1.);
  lpp.diameter_max_cutoff.push_back(-1.);
  lpp.density_particle[0] = 2500;
  lpp.number[0]           = 1000;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Calling normal distribution
  std::vector<std::shared_ptr<Distribution>> distribution_object_container;

  distribution_object_container.push_back(
    std::make_shared<LogNormalDistribution>(
      lpp.particle_average_diameter[0],
      lpp.particle_size_std[0],
      lpp.seed_for_distributions[0],
      lpp.diameter_min_cutoff[0],
      lpp.diameter_max_cutoff[0],
      lpp.distribution_weighting_type[0]));

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

  double             sum_dp                    = 0.;
  int                particle_number           = 0;
  const unsigned int total_number_of_particles = lpp.number[0];
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle, ++particle_number)
    {
      auto particle_properties = particle->get_properties();

      double dp = particle_properties[PropertiesIndex::dp];

      deallog << "Particle " << particle_number << " diameter is: " << dp
              << std::endl;

      sum_dp += dp;
    }
  const double mean_dp  = sum_dp / total_number_of_particles;
  double       variance = 0.;
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      auto   particle_properties = particle->get_properties();
      double dp                  = particle_properties[PropertiesIndex::dp];

      variance += std::pow(dp - mean_dp, 2);
    }
  variance                 = variance / total_number_of_particles;
  const double sigma_verif = std::sqrt(variance);

  deallog << "Distribution mean: " << mean_dp << std::endl;
  deallog << "Distribution standard distribution: " << sigma_verif << std::endl
          << std::endl;
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
           DistributionWeightingType::number_based>(0.005, 0.0005);

      // For the volume based test, the mu and sigma were choose to match
      // the same output as the number based one.
      test<3,
           DEM::DEMProperties::PropertiesIndex,
           DistributionWeightingType::volume_based>(0.005151505, 0.0005151505);
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
