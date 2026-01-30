// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later


#ifndef test_particles_functions_h
#define test_particles_functions_h


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

#include <dem/particle_interaction_outcomes.h>
#include <dem/particle_particle_contact_force.h>


// Tests (with common definitions)
#include <../tests/tests.h>


/**
 * @brief Return the particle iterator with the inserted particle, according to the position of the particle.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param triangulation Triangulation to access the information of the cells.
 * @param position Postion of the particle.
 * @param id Id of the particle.
 * @return Particle iterator with the inserted particle.
 */
template <int dim>
Particles::ParticleIterator<dim>
construct_particle_iterator(
  Particles::ParticleHandler<dim>           &particle_handler,
  parallel::distributed::Triangulation<dim> &triangulation,
  Point<3>                                  &position,
  int                                        id)
{
  typename Triangulation<dim>::active_cell_iterator cell =
    GridTools::find_active_cell_around_point(triangulation, position);
  Particles::Particle<dim> particle(position, position, id);
  return particle_handler.insert_particle(particle, cell);
}


/**
 * @brief the properties of the particle in the PropertiesIndex.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param pit Particle Iterator where the particle is inserted.
 * @param type Particle type.
 * @param diameter Particle diameter.
 * @param mass Particle mass.
 * @param v Particle velocity.
 * @param omega Particle angular velocity.
 */
template <int dim, typename PropertiesIndex>
void
set_particle_properties(Particles::ParticleIterator<dim> &pit,
                        int                               type,
                        double                            diameter,
                        double                            mass,
                        Tensor<1, dim>                   &v,
                        Tensor<1, dim>                   &omega)
{
  pit->get_properties()[PropertiesIndex::type]    = type;
  pit->get_properties()[PropertiesIndex::dp]      = diameter;
  pit->get_properties()[PropertiesIndex::mass]    = mass;
  pit->get_properties()[PropertiesIndex::v_x]     = v[0];
  pit->get_properties()[PropertiesIndex::omega_x] = omega[0];
  if (dim > 1)
    {
      pit->get_properties()[PropertiesIndex::v_y]     = v[1];
      pit->get_properties()[PropertiesIndex::omega_y] = omega[1];
    }
  if (dim > 2)
    {
      pit->get_properties()[PropertiesIndex::v_z]     = v[2];
      pit->get_properties()[PropertiesIndex::omega_z] = omega[2];
    }
}


/**
 * @brief Set default values for all lagrangian properties.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @param[in] particle_type_number Number of particle types.
 * @param[out] dem_parameters Simulation parameters
 */
template <int dim>
void
set_default_dem_parameters(const unsigned int        particle_type_number,
                           DEMSolverParameters<dim> &dem_parameters)
{
  Parameters::Lagrangian::LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;

  lpp.particle_type_number = particle_type_number;

  // Rolling resistance method
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant;

  // Particle parameters
  for (unsigned int i = 0; i < particle_type_number; ++i)
    {
      lpp.particle_average_diameter.push_back(0.1);
      lpp.custom_distribution_from_file.push_back(false);
      lpp.custom_distribution_filenames.emplace_back(" ");
      lpp.custom_probability_function_type.push_back(
        Parameters::Lagrangian::ProbabilityFunctionType::PDF);
      lpp.custom_distribution_interpolation.push_back(false);
      lpp.particle_custom_diameter.push_back(std::vector<double>{0.1, 0.2});
      lpp.particle_custom_probability.push_back(std::vector<double>{0.5, 0.5});
      lpp.distribution_weighting_type.push_back(
        Parameters::Lagrangian::DistributionWeightingType::number_based);
      lpp.seed_for_distributions.push_back(1);
      lpp.diameter_min_cutoff.push_back(-1.);
      lpp.diameter_max_cutoff.push_back(-1.);
      lpp.number.push_back(0);
      lpp.density_particle.push_back(1000);
      lpp.youngs_modulus_particle.push_back(1000000);
      lpp.poisson_ratio_particle.push_back(0.3);
      lpp.restitution_coefficient_particle.push_back(0.1);
      lpp.friction_coefficient_particle.push_back(0.1);
      lpp.rolling_viscous_damping_coefficient_particle.push_back(0.1);
      lpp.rolling_friction_coefficient_particle.push_back(0.1);
      lpp.surface_energy_particle.push_back(0.0);
      lpp.hamaker_constant_particle.push_back(4.e-19);
      lpp.thermal_conductivity_particle.push_back(1);
      lpp.specific_heat_particle.push_back(1000);
      lpp.microhardness_particle.push_back(1.e9);
      lpp.surface_slope_particle.push_back(0.1);
      lpp.surface_roughness_particle.push_back(1.e-9);
      lpp.thermal_accommodation_particle.push_back(0.7);
      lpp.real_youngs_modulus_particle.push_back(1.e9);
    }

  // Wall parameters
  lpp.youngs_modulus_wall          = 1000000;
  lpp.poisson_ratio_wall           = 0.3;
  lpp.restitution_coefficient_wall = 0.1;
  lpp.friction_coefficient_wall    = 0.1;
  lpp.rolling_friction_wall        = 0.1;
  lpp.rolling_viscous_damping_wall = 0.1;
  lpp.surface_energy_wall          = 0.0;
  lpp.hamaker_constant_wall        = 4.e-19;
  lpp.thermal_conductivity_wall    = 100;
  lpp.microhardness_wall           = 1.e9;
  lpp.surface_slope_wall           = 0.1;
  lpp.surface_roughness_wall       = 1.e-10;
  lpp.thermal_accommodation_wall   = 0.7;
  lpp.real_youngs_modulus_wall     = 1.e9;

  // Interstitial gas parameters
  lpp.thermal_conductivity_gas     = 0.01;
  lpp.specific_heat_gas            = 1000;
  lpp.dynamic_viscosity_gas        = 1.e-5;
  lpp.specific_heats_ratio_gas     = 1;
  lpp.molecular_mean_free_path_gas = 68.e-9;
}


/**
 * @brief Set the contact outcomes to 0 for each particle.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param contact_outcome Interaction outcomes.
 */
template <int dim, typename PropertiesIndex>
void
reinitialize_contact_outcomes(
  Particles::ParticleHandler<dim>              &particle_handler,
  ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Getting id of particle as local variable
      unsigned int particle_id = particle->get_local_index();

      // Reinitializing contact outcomes of particles in the system
      contact_outcome.force[particle_id][0] = 0;
      contact_outcome.force[particle_id][1] = 0;
      contact_outcome.force[particle_id][2] = 0;

      contact_outcome.torque[particle_id][0] = 0;
      contact_outcome.torque[particle_id][1] = 0;
      contact_outcome.torque[particle_id][2] = 0;

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::DEMMPProperties::PropertiesIndex>)
        {
          contact_outcome.heat_transfer_rate[particle_id] = 0;
        }
    }
}

#endif // test_particles_functions_h
