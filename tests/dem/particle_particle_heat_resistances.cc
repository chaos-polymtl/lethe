// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


/**
 * @brief This test checks the calculation of thermal resistances used to model conduction between two particles.
 *
 */


// Tests (with common definitions)
#include <../tests/dem/multiphysics_integrator.h>
#include <../tests/dem/test_particles_functions.h>
#include <../tests/dem/particle_particle_heat_transfer.h>


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
  DEMSolverParameters<dim> dem_parameters;
  Parameters::Lagrangian::LagrangianPhysicalProperties &lagrangian_prop =
    dem_parameters.lagrangian_physical_properties;
  Parameters::Lagrangian::ModelParameters &model_param =
    dem_parameters.model_parameters;

  const double specific_heat   = 840; // glass
  const Tensor<1, dim> g{{0, 0, 0}};
  const double         dt                                               = 0.001;
  const double         particle_diameter                                = 0.005;
  lagrangian_prop.particle_type_number                            = 1;
  lagrangian_prop.youngs_modulus_particle[0]                      = 5000000;
  lagrangian_prop.poisson_ratio_particle[0]                       = 0.22;
  lagrangian_prop.restitution_coefficient_particle[0]             = 0.8;
  lagrangian_prop.friction_coefficient_particle[0]                = 1;
  lagrangian_prop.rolling_friction_coefficient_particle[0]        = 0.02;
  lagrangian_prop.density_particle[0]                             = 2521;
  model_param.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  const double neighborhood_threshold = std::pow(1.3 * particle_diameter, 2);

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
  Point<3> position_1 = {0.4, 0, 0};
  const unsigned int      id_1       = 0;
  Point<3> position_2 = {0.40499, 0, 0};
  const unsigned int      id_2       = 1;

  // Constructing particle iterators from particle positions (inserting
  // particles)
  Particles::ParticleIterator<dim> pit_1 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position_1, id_1);
  Particles::ParticleIterator<dim> pit_2 = construct_particle_iterator<dim>(
    particle_handler, triangulation, position_2, id_2);

  // Setting particle properties
  Tensor<1, dim> v_1{{0, 0, 0}};
  Tensor<1, dim> omega_1{{0, 0, 0}};
  Tensor<1, dim> v_2{{0, 0, 0}};
  Tensor<1, dim> omega_2{{0, 0, 0}};
  const double         mass = 1;
  const unsigned int            type = 0;
  const double T_initial_1 = 300;
  const double T_initial_2 = 100;

  set_particle_properties<dim, PropertiesIndex>(
    pit_1, type, particle_diameter, mass, v_1, omega_1);
  set_particle_properties<dim, PropertiesIndex>(
    pit_2, type, particle_diameter, mass, v_2, omega_2);

  pit_1->get_properties()[PropertiesIndex::T]             = T_initial_1;
  pit_1->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  pit_2->get_properties()[PropertiesIndex::T]             = T_initial_2;
  pit_2->get_properties()[PropertiesIndex::specific_heat] = specific_heat;

  // Initialize variables
  std::vector<Tensor<1, 3>> torque;
  std::vector<Tensor<1, 3>> force;
  std::vector<double>       MOI;
  std::vector<double> heat_transfer;
  std::vector<double> heat_source;

  particle_handler.sort_particles_into_subdomains_and_cells();
  force.resize(particle_handler.get_max_local_particle_index());
  torque.resize(force.size());
  MOI.resize(force.size());
  for (auto &moi_val : MOI)
    moi_val = 1;
  heat_transfer.resize(force.size());
  heat_source.resize(force.size());

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
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
    nonlinear_force_object(dem_parameters);
  nonlinear_force_object.calculate_particle_particle_contact_force(
    contact_manager.get_local_adjacent_particles(),
    contact_manager.get_ghost_adjacent_particles(),
    contact_manager.get_local_local_periodic_adjacent_particles(),
    contact_manager.get_local_ghost_periodic_adjacent_particles(),
    contact_manager.get_ghost_local_periodic_adjacent_particles(),
    dt,
    torque,
    force);

  auto particle_one = particle_handler.begin();
  auto particle_two = std::next(particle_one);
  auto particle_one_properties = particle_one->get_properties();
  auto particle_two_properties = particle_two->get_properties();
  Point<3> particle_one_location = particle_one->get_location();
  Point<3> particle_two_location = particle_two->get_location();

  const double radius_one = particle_diameter / 2 ;
  const double radius_two = particle_diameter / 2 ;
  const double effective_radius = (radius_one*radius_two) / (radius_one+radius_two);
  const double effective_youngs_modulus = 1/2 * lagrangian_prop.youngs_modulus_particle.at(0);
  const double effective_real_youngs_modulus = 65.e9;
  const double effective_surface_roughness = 25e-9;
  const double effective_surface_slope = 0.078;
  const double effective_microhardness = 9e9;
  const double thermal_conductivity_one = 1;
  const double thermal_conductivity_two = 1;
  const double thermal_conductivity_gas = 0.027;
  const double gas_parameter_m = 3.66e-7;
  const double normal_overlap =
  0.5 * (particle_one_properties[PropertiesIndex::dp] +
         particle_two_properties[PropertiesIndex::dp]) -
  particle_one_location.distance(particle_two_location);
  const double normal_force_norm = force[particle_one->get_id()].norm();

  const double effective_particle_conductivity =
    (thermal_conductivity_one * thermal_conductivity_two) /
    (thermal_conductivity_one + thermal_conductivity_two);

  // Calculation of contact radius
  const double contact_radius =
    calculate_corrected_contact_radius(effective_radius,
                                       effective_youngs_modulus,
                                       effective_real_youngs_modulus,
                                       normal_force_norm);
  const double maximum_pressure =
    (2 * effective_real_youngs_modulus * normal_overlap) /
    (M_PI * contact_radius);

  // Calculation of each thermal resistance
  const double R_L =
    calculate_macrocontact_spreading_resistance(effective_particle_conductivity,
                                                contact_radius);

  const double R_s =
    calculate_microcontact_spreading_resistance(effective_surface_slope,
                                                effective_surface_roughness,
                                                effective_microhardness,
                                                contact_radius,
                                                effective_particle_conductivity,
                                                maximum_pressure);

  const double R_c =
    calculate_macrogap_solid_layers_resistance(radius_one,
                                               radius_two,
                                               thermal_conductivity_one,
                                               thermal_conductivity_two,
                                               contact_radius);

  const double R_g =
    calculate_microgap_interstitial_gas_resistance(effective_surface_roughness,
                                                   contact_radius,
                                                   gas_parameter_m,
                                                   thermal_conductivity_gas,
                                                   maximum_pressure,
                                                   effective_microhardness);

  const double R_G =
    calculate_macrogap_interstitial_gas_resistance(effective_radius,
                                                   thermal_conductivity_gas,
                                                   contact_radius,
                                                   gas_parameter_m);

    




 
  // Output
  deallog << contact_radius << " " <<  R_L << " " << R_s << " " << R_c << " " << R_g << " " << R_G << std::endl;
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
