// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
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
 * @param particle_diameter Particle diameter.
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
  Parameters::Lagrangian::LagrangianPhysicalProperties &properties =
    dem_parameters.lagrangian_physical_properties;

  properties.particle_type_number = particle_type_number;

  // Rolling resistance method
  dem_parameters.model_parameters.rolling_resistance_method =
    Parameters::Lagrangian::RollingResistanceMethod::constant_resistance;

  // Particle parameters
  for (unsigned int i = 0; i < particle_type_number; ++i)
    {
      properties.density_particle[i]                             = 1000;
      properties.youngs_modulus_particle[i]                      = 1000000;
      properties.poisson_ratio_particle[i]                       = 0.3;
      properties.restitution_coefficient_particle[i]             = 0.1;
      properties.friction_coefficient_particle[i]                = 0.1;
      properties.rolling_friction_coefficient_particle[i]        = 0.1;
      properties.rolling_viscous_damping_coefficient_particle[i] = 0.1;
      properties.surface_energy_particle[i]                      = 0.0;
      properties.hamaker_constant_particle[i]                    = 4.e-19;
      properties.thermal_conductivity_particle[i]                = 1;
      properties.specific_heat_particle[i]                       = 1000;
      properties.microhardness_particle[i]                       = 1.e9;
      properties.surface_slope_particle[i]                       = 0.1;
      properties.surface_roughness_particle[i]                   = 1.e-9;
      properties.thermal_accommodation_particle[i]               = 0.7;
      properties.real_youngs_modulus_particle[i]                 = 1.e9;
    }

  // Wall parameters
  properties.youngs_modulus_wall          = 1000000;
  properties.poisson_ratio_wall           = 0.3;
  properties.restitution_coefficient_wall = 0.1;
  properties.friction_coefficient_wall    = 0.1;
  properties.rolling_friction_wall        = 0.1;
  properties.rolling_viscous_damping_wall = 0.1;
  properties.surface_energy_wall          = 0.0;
  properties.hamaker_constant_wall        = 4.e-19;
  properties.thermal_conductivity_wall    = 100;
  properties.microhardness_wall           = 1.e9;
  properties.surface_slope_wall           = 0.1;
  properties.surface_roughness_wall       = 1.e-10;
  properties.thermal_accommodation_wall   = 0.7;
  properties.real_youngs_modulus_wall     = 1.e9;

  // Interstitial gas parameters
  properties.thermal_conductivity_gas     = 0.01;
  properties.specific_heat_gas            = 1000;
  properties.dynamic_viscosity_gas        = 1.e-5;
  properties.specific_heats_ratio_gas     = 1;
  properties.molecular_mean_free_path_gas = 68.e-9;
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
      unsigned int particle_id = particle->get_id();

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
