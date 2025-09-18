// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


/**
 * @brief This test
 *
 */

#include <../tests/dem/test_particles_functions.h>
#include <dem/adaptive_sparse_contacts.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_particle_contact_force.h>

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  const int                                 hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  const int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining simulation general parameters
  DEMSolverParameters<dim>                              dem_parameters;
  Parameters::Lagrangian::LagrangianPhysicalProperties &properties =
    dem_parameters.lagrangian_physical_properties;
  Parameters::Lagrangian::ModelParameters<dim> &model_param =
    dem_parameters.model_parameters;

  const Tensor<1, dim> g{{0, 0, 0}};
  const double         dt                                    = 0.1;
  const double         particle_diameter                     = 0.01;
  const double         specific_heat                         = 840;
  properties.particle_type_number                            = 1;
  properties.youngs_modulus_particle[0]                      = 65e9;
  properties.poisson_ratio_particle[0]                       = 0.22;
  properties.restitution_coefficient_particle[0]             = 0.8;
  properties.friction_coefficient_particle[0]                = 1;
  properties.rolling_friction_coefficient_particle[0]        = 0.02;
  properties.density_particle[0]                             = 2521;
  properties.rolling_viscous_damping_coefficient_particle[0] = 0.;
  properties.surface_energy_particle[0]                      = 0.;
  properties.hamaker_constant_particle[0]                    = 0.;
  model_param.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  // Defining parameters for thermal DEM
  properties.surface_roughness_particle[0]     = 25e-9;
  properties.surface_slope_particle[0]         = 0.078;
  properties.microhardness_particle[0]         = 9e9;
  properties.thermal_conductivity_particle[0]  = 1;
  properties.thermal_conductivity_gas          = 0.027;
  properties.dynamic_viscosity_gas             = 1.85e-5;
  properties.specific_heat_gas                 = 1006;
  properties.specific_heats_ratio_gas          = 1;
  properties.molecular_mean_free_path_gas      = 68e-9;
  properties.thermal_accommodation_particle[0] = 0.7;
  properties.real_youngs_modulus_particle[0]   = 65e9;

  // Defining particle handler
  Particles::ParticleHandler<dim> particle_handler(
    triangulation, mapping, PropertiesIndex::n_properties);

  // Creating containers manager for finding cell neighbor and also broad and
  // fine particle-particle search objects
  DEMContactManager<dim, PropertiesIndex> contact_manager;

  // Finding cell neighbors
  typename dem_data_structures<dim>::periodic_boundaries_cells_info
    dummy_pbc_info;
  contact_manager.execute_cell_neighbors_search(triangulation, dummy_pbc_info);

  // Inserting two particles in contact
  Point<3>           position_1 = {0, 0, 0};
  const unsigned int id_1       = 0;
  Point<3>           position_2 = {0.00999, 0, 0};
  const unsigned int id_2       = 1;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit_1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position_1, id_1);
  Particles::ParticleIterator<dim> pit_2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position_2, id_2);

  // Setting particle properties
  Tensor<1, dim>     v_1{{0, 0, 0}};
  Tensor<1, dim>     omega_1{{0, 0, 0}};
  Tensor<1, dim>     v_2{{0, 0, 0}};
  Tensor<1, dim>     omega_2{{0, 0, 0}};
  const double       mass        = 1;
  const unsigned int type        = 0;
  const double       T_initial_1 = 600;
  const double       T_initial_2 = 100;

  set_particle_properties<dim, PropertiesIndex>(
    pit_1, type, particle_diameter, mass, v_1, omega_1);
  set_particle_properties<dim, PropertiesIndex>(
    pit_2, type, particle_diameter, mass, v_2, omega_2);

  pit_1->get_properties()[PropertiesIndex::T]             = T_initial_1;
  pit_1->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  pit_2->get_properties()[PropertiesIndex::T]             = T_initial_2;
  pit_2->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  // Initializing variables
  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  std::vector<double>                          MOI;

  particle_handler.sort_particles_into_subdomains_and_cells();
  const unsigned int number_of_particles =
    particle_handler.get_max_local_particle_index();
  contact_outcome.resize_interaction_containers(number_of_particles);
  MOI.resize(number_of_particles);
  for (auto &moi_val : MOI)
    moi_val = 1;

  contact_manager.update_local_particles_in_cells(particle_handler);

  // Dummy Adaptive sparse contacts object and particle-particle broad search
  AdaptiveSparseContacts<dim, PropertiesIndex> dummy_adaptive_sparse_contacts;
  contact_manager.execute_particle_particle_broad_search(
    particle_handler, dummy_adaptive_sparse_contacts);

  // Calling fine search
  contact_manager.execute_particle_particle_fine_search(neighborhood_threshold);

  // Calculating and applying contact force and heat transfer rate
  ParticleParticleContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleParticleContactForceModel::
      hertz_mindlin_limit_overlap,
    Parameters::Lagrangian::RollingResistanceMethod::constant>
    nonlinear_force_object(dem_parameters);
  nonlinear_force_object.calculate_particle_particle_contact(
    contact_manager.get_local_adjacent_particles(),
    contact_manager.get_ghost_adjacent_particles(),
    contact_manager.get_local_local_periodic_adjacent_particles(),
    contact_manager.get_local_ghost_periodic_adjacent_particles(),
    contact_manager.get_ghost_local_periodic_adjacent_particles(),
    dt,
    contact_outcome);

  // Output
  auto particle_one = particle_handler.begin();
  deallog << "The heat transfer applied to particle one is "
          << contact_outcome.heat_transfer_rate[particle_one->get_id()]
          << " J/s." << std::endl;
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
