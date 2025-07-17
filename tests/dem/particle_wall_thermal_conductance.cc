// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the calculation of particle-wall
 * thermal conductance is checked.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>

#include <dem/contact_info.h>
#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_heat_transfer.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>

// Tests (with common definitions)
#include <../tests/dem/test_particles_functions.h>

#include <../tests/tests.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  const int                                 hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  const double grid_radius       = 0.5 * GridTools::diameter(tr);
  const int    refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);


  // Defining general simulation parameters
  DEMSolverParameters<dim> dem_parameters;
  set_default_dem_parameters(1, dem_parameters);
  auto &properties = dem_parameters.lagrangian_physical_properties;
  const Tensor<1, dim> g{{0, 0, -9.81}};
  const double         dt                                    = 1.e-5;
  const double         particle_diameter                     = 0.005;
  const unsigned int   rotating_wall_maximum_number          = 6;
  const double         poisson_ratio                         = 0.3;
  const double         youngs_modulus                        = 5.e6;
  properties.particle_type_number                            = 1;
  properties.youngs_modulus_particle[0]                      = youngs_modulus;
  properties.youngs_modulus_wall                             = youngs_modulus;
  properties.poisson_ratio_particle[0]                       = poisson_ratio;
  properties.poisson_ratio_wall                              = poisson_ratio;
  properties.restitution_coefficient_particle[0]             = 0.5;
  properties.restitution_coefficient_wall                    = 0.5;
  properties.friction_coefficient_particle[0]                = 0.5;
  properties.friction_coefficient_wall                       = 0.5;
  properties.rolling_friction_coefficient_particle[0]        = 0.1;
  properties.rolling_friction_wall                           = 0.1;
  properties.rolling_viscous_damping_coefficient_particle[0] = 0.1;
  properties.rolling_viscous_damping_wall                    = 0.1;
  properties.density_particle[0]                             = 2500;
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  // Defining parameters for thermal DEM
  const double equivalent_surface_roughness  = 1e-9;
  const double equivalent_surface_slope      = 0.08;
  const double effective_microhardness       = 9e9;
  const double thermal_conductivity_particle = 3000;
  const double thermal_conductivity_wall     = 300;
  const double thermal_conductivity_gas      = 0.2;
  const double dynamic_viscosity_gas         = 9.e-6;
  const double specific_heat_gas             = 10000;
  const double specific_heats_ratio_gas      = 1.4;
  const double molecular_mean_free_path_gas  = 68e-9;
  const double thermal_accommodation         = 0.7;

  // Initializing motion of boundaries
  Tensor<1, dim> translational_and_rotational_velocity;
  for (unsigned int d = 0; d < dim; ++d)
    {
      translational_and_rotational_velocity[d] = 0;
    }
  for (unsigned int counter = 0; counter < rotating_wall_maximum_number;
       ++counter)
    {
      dem_parameters.boundary_conditions.boundary_rotational_speed.insert(
        {counter, 0});
      dem_parameters.boundary_conditions.boundary_translational_velocity.insert(
        {counter, translational_and_rotational_velocity});
      dem_parameters.boundary_conditions.boundary_rotational_vector.insert(
        {counter, translational_and_rotational_velocity});
    }

  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Inserting one particle in contact with a wall
  Point<dim>                       position_1 = {-0.998, 0, 0};
  const unsigned int               id_1       = 0;
  Particles::ParticleIterator<dim> pit_1 =
    construct_particle_iterator<dim>(particle_handler, tr, position_1, id_1);

  // Setting particle properties
  Tensor<1, dim>     v_1{{0.01, 0, 0}};
  Tensor<1, dim>     omega_1{{0, 0, 0}};
  const double       mass = 1;
  const unsigned int type = 0;
  set_particle_properties<dim, PropertiesIndex>(
    pit_1, type, particle_diameter, mass, v_1, omega_1);

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

  // Finding boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::vector<unsigned int>     outlet_boundaries;
  boundary_cells_object.build(
    tr,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Calling broad search
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    particle_wall_contact_list;
  find_particle_wall_contact_pairs<dim>(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);
  find_particle_wall_contact_pairs<dim>(
    boundary_cells_object.get_boundary_cells_information(),
    particle_handler,
    particle_wall_contact_list);

  // Calling fine search
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_pairs_in_contact;
  particle_wall_fine_search<dim>(particle_wall_contact_list,
                                 particle_wall_pairs_in_contact);

  // Calling non-linear force
  ParticleWallContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleWallContactForceModel::nonlinear,
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance>
                           nonlinear_force_object(dem_parameters);
  bool                     floating_wall = false;
  OngoingCollisionLog<dim> ongoing_collision_log;
  CollisionEventLog<dim>   collision_event_log;
  double                   time;
  nonlinear_force_object.calculate_particle_wall_contact(
    particle_wall_pairs_in_contact,
    dt,
    time,
    particle_handler,
    floating_wall,
    contact_outcome,
    ongoing_collision_log,
    collision_event_log);

  // Calculating  parameters for thermal DEM
  auto     particle          = particle_handler.begin();
  Point<3> particle_location = particle->get_location();
  auto    &pairs_in_contact_content =
    particle_wall_pairs_in_contact.at(particle->get_id());
  auto &contact_info = pairs_in_contact_content.begin()->second;

  const double radius_particle = particle_diameter * 0.5;
  const double effective_youngs_modulus =
    youngs_modulus / (2.0 * (1. - poisson_ratio * poisson_ratio));
  const Tensor<1, 3> point_to_particle_vector =
    particle_location - contact_info.point_on_boundary;
  const Tensor<1, 3> normal_vector = contact_info.normal_vector;
  const Tensor<1, 3> projected_vector =
    ((point_to_particle_vector * normal_vector) /
     (normal_vector.norm_square())) *
    normal_vector;
  const double normal_overlap = radius_particle - (projected_vector.norm());
  const double normal_force_norm =
    contact_outcome.force[particle->get_id()].norm();
  const double prandtl_gas =
    dynamic_viscosity_gas * specific_heat_gas / thermal_conductivity_gas;
  const double gas_parameter_m =
    2. * (2. - thermal_accommodation) / thermal_accommodation *
    (2. * specific_heats_ratio_gas) / (1. + specific_heats_ratio_gas) *
    molecular_mean_free_path_gas / prandtl_gas;

  double thermal_conductance;
  calculate_contact_thermal_conductance<ContactType::particle_floating_mesh>(
    radius_particle,
    0, // unused
    effective_youngs_modulus,
    effective_youngs_modulus,
    equivalent_surface_roughness,
    equivalent_surface_slope,
    effective_microhardness,
    thermal_conductivity_particle,
    thermal_conductivity_wall,
    thermal_conductivity_gas,
    gas_parameter_m,
    normal_overlap,
    normal_force_norm,
    thermal_conductance);

  // Output
  deallog << "The contact thermal conductance is: " << thermal_conductance
          << std::endl;
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
