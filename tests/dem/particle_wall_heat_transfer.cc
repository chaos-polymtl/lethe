// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief In this test, the heat transfer between particles and the walls of
 * the grid is checked. The DEM boundary conditions are parsed from a parameter
 * string. Particles 0 and 2 are in contact with an isothermal wall whose
 * temperature depends on time and space, while particle 1 is in contact with
 * an adiabatic wall.
 *
 * The contact configuration and the physical properties of particles 0 and 2
 * are the same as in the particle_wall_thermal_conductance test, so their
 * thermal conductance H is the one obtained in that test. The temperature of
 * the wall differs from the temperature of the particles by -1 K or 1 K for
 * particle 0, and by 1 K or 3 K for particle 2, which touches a hotter part of
 * the wall. Hence, the heat transfer rates must be -H, H or 3H. The heat
 * transfer rate of particle 1 must be zero.
 */

// Deal.II
#include <deal.II/base/parameter_handler.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters_lagrangian.h>

#include <dem/dem_solver_parameters.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_fine_search.h>

// Tests (with common definitions)
#include <../tests/dem/test_particles_functions.h>

#include <../tests/tests.h>

#include <algorithm>
#include <map>
#include <set>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
test()
{
  // Creating the mesh and refinement. The boundaries are colorized: boundary 0
  // is the wall at x = -1 and boundary 2 is the wall at y = -1.
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);
  const int                                 hyper_cube_length = 1;
  GridGenerator::hyper_cube(tr,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  const int refinement_number = 2;
  tr.refine_global(refinement_number);
  MappingQ<dim> mapping(1);

  // Defining general simulation parameters
  DEMSolverParameters<dim> dem_parameters;
  set_default_dem_parameters(1, dem_parameters);
  auto &properties = dem_parameters.lagrangian_physical_properties;

  const double dt                = 1.e-5;
  const double particle_diameter = 0.005;
  const double poisson_ratio     = 0.3;
  const double youngs_modulus    = 5.e6;

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

  // Defining parameters for thermal DEM. They are chosen such that the
  // effective particle-wall properties are the same as in the
  // particle_wall_thermal_conductance test.
  properties.real_youngs_modulus_particle[0]   = youngs_modulus;
  properties.real_youngs_modulus_wall          = youngs_modulus;
  properties.surface_roughness_particle[0]     = 1.e-9;
  properties.surface_roughness_wall            = 0.;
  properties.surface_slope_particle[0]         = 0.08;
  properties.surface_slope_wall                = 0.;
  properties.microhardness_particle[0]         = 9.e9;
  properties.microhardness_wall                = 9.e9;
  properties.thermal_conductivity_particle[0]  = 3000;
  properties.thermal_conductivity_wall         = 300;
  properties.thermal_accommodation_particle[0] = 0.7;
  properties.thermal_accommodation_wall        = 0.7;
  properties.thermal_conductivity_gas          = 0.2;
  properties.dynamic_viscosity_gas             = 9.e-6;
  properties.specific_heat_gas                 = 10000;
  properties.specific_heats_ratio_gas          = 1.4;
  properties.molecular_mean_free_path_gas      = 68e-9;

  // Parsing the DEM boundary conditions. The wall at x = -1 is isothermal: its
  // temperature is 19 before t = 0.5 and 21 after, and it is 2 degrees higher
  // where y > 0.1. The wall at y = -1 is adiabatic.
  ParameterHandler prm;
  dem_parameters.boundary_conditions.declare_parameters(prm);
  prm.parse_input_from_string(R"(
subsection DEM boundary conditions
  set number of boundary conditions = 2
  subsection boundary condition 0
    set boundary id           = 0
    set type                  = fixed_wall
    set thermal boundary type = isothermal
    subsection wall temperature
      set Function expression = if(t < 0.5, 19, 21) + if(y > 0.1, 2, 0)
    end
  end
  subsection boundary condition 1
    set boundary id           = 2
    set type                  = fixed_wall
    set thermal boundary type = adiabatic
  end
end
)");
  dem_parameters.boundary_conditions.parse_parameters(prm);

  Particles::ParticleHandler<dim> particle_handler(
    tr, mapping, PropertiesIndex::n_properties);

  // Inserting particles 0 and 2 in contact with the isothermal wall, at y = 0
  // and y = 0.25, and particle 1 in contact with the adiabatic wall. The
  // particles move away from their wall in the normal direction, so that there
  // is no tangential force.
  const double       mass = 1;
  const unsigned int type = 0;
  Tensor<1, dim>     omega{{0, 0, 0}};

  Point<dim>                       position_0 = {-0.998, 0, 0};
  Particles::ParticleIterator<dim> pit_0 =
    construct_particle_iterator<dim>(particle_handler, tr, position_0, 0);
  Tensor<1, dim> v_0{{0.01, 0, 0}};
  set_particle_properties<dim, PropertiesIndex>(
    pit_0, type, particle_diameter, mass, v_0, omega);
  pit_0->get_properties()[PropertiesIndex::T] = 20;

  Point<dim>                       position_1 = {0.25, -0.998, 0.25};
  Particles::ParticleIterator<dim> pit_1 =
    construct_particle_iterator<dim>(particle_handler, tr, position_1, 1);
  Tensor<1, dim> v_1{{0, 0.01, 0}};
  set_particle_properties<dim, PropertiesIndex>(
    pit_1, type, particle_diameter, mass, v_1, omega);
  pit_1->get_properties()[PropertiesIndex::T] = 0;

  Point<dim>                       position_2 = {-0.998, 0.25, 0.25};
  Particles::ParticleIterator<dim> pit_2 =
    construct_particle_iterator<dim>(particle_handler, tr, position_2, 2);
  Tensor<1, dim> v_2{{0.01, 0, 0}};
  set_particle_properties<dim, PropertiesIndex>(
    pit_2, type, particle_diameter, mass, v_2, omega);
  pit_2->get_properties()[PropertiesIndex::T] = 20;

  // Initializing variables
  ParticleInteractionOutcomes<PropertiesIndex> contact_outcome;
  particle_handler.sort_particles_into_subdomains_and_cells();
  contact_outcome.resize_interaction_containers(
    particle_handler.get_max_local_particle_index());

  // Finding boundary cells
  BoundaryCellsInformation<dim> boundary_cells_object;
  std::set<types::boundary_id>  outlet_boundaries;
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

  // Calling fine search
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_pairs_in_contact;
  particle_wall_fine_search<dim>(particle_wall_contact_list,
                                 particle_wall_pairs_in_contact);

  // Non-linear particle-wall contact force object
  ParticleWallContactForce<
    dim,
    PropertiesIndex,
    Parameters::Lagrangian::ParticleWallContactForceModel::nonlinear,
    Parameters::Lagrangian::RollingResistanceMethod::constant>
    nonlinear_force_object(dem_parameters);

  // The particles do not move, so the contact configuration is the same at
  // both times and only the temperature of the isothermal wall changes.
  for (const double time : {0., 1.})
    {
      std::ranges::fill(contact_outcome.heat_transfer_rate, 0.);

      nonlinear_force_object.update_boundary_temperature(time);
      nonlinear_force_object.calculate_particle_wall_contact(
        particle_wall_pairs_in_contact, dt, contact_outcome);

      // Gather the heat transfer rate of each particle, sorted by id
      std::map<types::particle_index, double> heat_transfer_rates;
      for (auto particle = particle_handler.begin();
           particle != particle_handler.end();
           ++particle)
        heat_transfer_rates[particle->get_id()] =
          contact_outcome.heat_transfer_rate[particle->get_local_index()];

      // Output
      deallog << "Time: " << time << std::endl;
      deallog << "Heat transfer rate of particle 0 (isothermal wall, y = 0): "
              << heat_transfer_rates.at(0) << std::endl;
      deallog << "Heat transfer rate of particle 1 (adiabatic wall): "
              << heat_transfer_rates.at(1) << std::endl;
      deallog
        << "Heat transfer rate of particle 2 (isothermal wall, y = 0.25): "
        << heat_transfer_rates.at(2) << std::endl;
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
