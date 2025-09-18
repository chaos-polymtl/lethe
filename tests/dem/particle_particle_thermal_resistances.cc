// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


/**
 * @brief This test checks the calculation of the thermal resistance network used to model conduction between two particles.
 *
 */


#include <../tests/dem/test_particles_functions.h>
#include <dem/dem_contact_manager.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_particle_contact_force.h>

#include <deal.II/base/logstream.h>

#include <../tests/tests.h>

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

  set_default_dem_parameters(1, dem_parameters);

  const Tensor<1, dim> g{{0, 0, 0}};
  const double         dt                             = 0.001;
  const double         particle_diameter              = 0.01;
  const double         youngs_modulus                 = 5000000;
  const double         poisson_ratio                  = 0.22;
  const double         specific_heat                  = 840;
  properties.particle_type_number                     = 1;
  properties.youngs_modulus_particle[0]               = youngs_modulus;
  properties.poisson_ratio_particle[0]                = poisson_ratio;
  properties.restitution_coefficient_particle[0]      = 0.8;
  properties.friction_coefficient_particle[0]         = 1;
  properties.rolling_friction_coefficient_particle[0] = 0.02;
  properties.density_particle[0]                      = 2521;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

  // Defining parameters for thermal DEM
  const double real_youngs_modulus          = 65.e9;
  const double equivalent_surface_roughness = 25e-9;
  const double equivalent_surface_slope     = 0.078;
  const double effective_microhardness      = 9e9;
  const double thermal_conductivity_one     = 1;
  const double thermal_conductivity_two     = 1;
  const double thermal_conductivity_gas     = 0.027;
  const double dynamic_viscosity_gas        = 1.85e-5;
  const double specific_heat_gas            = 1006;
  const double specific_heats_ratio_gas     = 1;
  const double molecular_mean_free_path_gas = 68e-9;
  const double thermal_accommodation        = 0.7;

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
  Point<3>           position_1 = {0., 0, 0};
  const unsigned int id_1       = 0;
  Point<3>           position_2 = {0.0099, 0, 0};
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
  const double       T_initial_1 = 300;
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

  // Calling non linear force
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

  // Calculating additional parameters for thermal DEM
  auto     particle_one          = particle_handler.begin();
  auto     particle_two          = std::next(particle_one);
  Point<3> particle_one_location = particle_one->get_location();
  Point<3> particle_two_location = particle_two->get_location();

  const double radius_one = particle_diameter * 0.5;
  const double radius_two = particle_diameter * 0.5;
  const double effective_youngs_modulus =
    youngs_modulus / (2.0 * (1. - poisson_ratio * poisson_ratio));
  const double effective_real_youngs_modulus =
    real_youngs_modulus / (2.0 * (1. - poisson_ratio * poisson_ratio));
  const double harmonic_particle_conductivity =
    harmonic_mean(thermal_conductivity_one, thermal_conductivity_two);
  const double harmonic_radius = harmonic_mean(radius_one, radius_two);
  const double normal_overlap =
    (radius_one + radius_two) -
    particle_one_location.distance(particle_two_location);
  const double normal_force_norm =
    contact_outcome.force[particle_one->get_id()].norm();
  const double prandtl_gas =
    dynamic_viscosity_gas * specific_heat_gas / thermal_conductivity_gas;
  const double gas_parameter_m =
    2. * (2. - thermal_accommodation) / thermal_accommodation *
    (2. * specific_heats_ratio_gas) / (1. + specific_heats_ratio_gas) *
    molecular_mean_free_path_gas / prandtl_gas;

  // Calculating the contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(harmonic_radius * 0.5,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);
  const double contact_radius_squared = contact_radius * contact_radius;

  const double corrected_normal_overlap =
    normal_overlap *
    pow(effective_youngs_modulus / effective_real_youngs_modulus, 2.0 / 3.0);
  const double maximum_pressure =
    (2.0 * effective_real_youngs_modulus * corrected_normal_overlap) /
    (M_PI * contact_radius);

  // Calculating each thermal resistance
  const double resistance_macrocontact =
    calculate_macrocontact_resistance(harmonic_particle_conductivity,
                                      contact_radius);

  const double resistance_microcontact =
    calculate_microcontact_resistance(equivalent_surface_slope,
                                      equivalent_surface_roughness,
                                      effective_microhardness,
                                      contact_radius_squared,
                                      harmonic_particle_conductivity,
                                      maximum_pressure);

  const double resistance_solid_macrogap =
    calculate_solid_macrogap_resistance(radius_one,
                                        thermal_conductivity_one,
                                        contact_radius_squared) +
    calculate_solid_macrogap_resistance(radius_two,
                                        thermal_conductivity_two,
                                        contact_radius_squared);

  const double resistance_gas_microgap =
    calculate_interstitial_gas_microgap_resistance(equivalent_surface_roughness,
                                                   contact_radius_squared,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  const double resistance_gas_macrogap =
    calculate_interstitial_gas_macrogap_resistance(harmonic_radius,
                                                   thermal_conductivity_gas,
                                                   contact_radius_squared,
                                                   gas_parameter_m);

  // Calculating total conductance and total resistance
  const double total_conductance =
    1. / (resistance_macrocontact +
          1. / (1. / resistance_microcontact + 1. / resistance_gas_microgap)) +
    1. / (resistance_solid_macrogap + resistance_gas_macrogap);
  const double total_resistance = 1. / total_conductance;

  // Output
  deallog << " contact radius : " << contact_radius << std::endl
          << " resistance macrocontact : " << resistance_macrocontact
          << std::endl
          << " resistance microcontact: " << resistance_microcontact
          << std::endl
          << " resistance solid macrogap : " << resistance_solid_macrogap
          << std::endl
          << " resistance gas microgap : " << resistance_gas_microgap
          << std::endl
          << " resistance gas macrogap : " << resistance_gas_macrogap
          << std::endl
          << " total resistance : " << total_resistance << std::endl
          << " total conductance : " << total_conductance << std::endl;
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
