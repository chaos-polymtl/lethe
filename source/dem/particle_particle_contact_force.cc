// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_lagrangian.h>

#include <dem/particle_particle_contact_force.h>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
ParticleParticleContactForce<dim,
                             PropertiesIndex,
                             contact_model,
                             rolling_friction_model>::
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
  , f_coefficient_epsd(dem_parameters.model_parameters.f_coefficient_epsd)
{
  set_effective_properties(dem_parameters);
  if constexpr (std::is_same_v<PropertiesIndex,
                               DEM::DEMMPProperties::PropertiesIndex>)
    {
      set_multiphysic_properties(dem_parameters);
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticleParticleContactForce<dim,
                             PropertiesIndex,
                             contact_model,
                             rolling_friction_model>::
  calculate_particle_particle_contact_force_and_heat_transfer_rate(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_local_periodic_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_ghost_periodic_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
                              &ghost_local_periodic_adjacent_particles,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    std::vector<double>       &heat_transfer_rate)
{
  // Calculating the contact forces and heat transfer rates for local-local
  // adjacent particles.
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::local_particle_particle>(
        adjacent_particles_list, torque, force, heat_transfer_rate, dt);
    }

  // Calculating the contact forces and heat transfer rates for local-ghost
  // adjacent particles.
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::ghost_particle_particle>(
        adjacent_particles_list, torque, force, heat_transfer_rate, dt);
    }

  // Calculating the contact forces and heat transfer rates for local-local
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::local_periodic_particle_particle>(
        periodic_adjacent_particles_list,
        torque,
        force,
        heat_transfer_rate,
        dt);
    }

  // Calculating the contact forces and heat transfer rates for local-ghost
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_ghost_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_periodic_particle_particle>(
        periodic_adjacent_particles_list,
        torque,
        force,
        heat_transfer_rate,
        dt);
    }

  // Calculating the contact forces and heat transfer rates for ghost-local
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       ghost_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_local_periodic_particle_particle>(
        periodic_adjacent_particles_list,
        torque,
        force,
        heat_transfer_rate,
        dt);
    }
}

template <int dim,
          typename PropertiesIndex,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticleParticleContactForce<dim,
                             PropertiesIndex,
                             contact_model,
                             rolling_friction_model>::
  calculate_particle_particle_contact_force(
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_local_periodic_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
      &local_ghost_periodic_adjacent_particles,
    typename dem_data_structures<dim>::adjacent_particle_pairs
                              &ghost_local_periodic_adjacent_particles,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force)
{
  std::vector<double> heat_transfer_rate;

  // Calling calculate_particle_particle_contact_force_and_heat_transfer_rate
  // with empty heat_transfer_rate
  calculate_particle_particle_contact_force_and_heat_transfer_rate(
    local_adjacent_particles,
    ghost_adjacent_particles,
    local_local_periodic_adjacent_particles,
    local_ghost_periodic_adjacent_particles,
    ghost_local_periodic_adjacent_particles,
    dt,
    torque,
    force,
    heat_transfer_rate);
}

// dem
// No resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;

// cfd_dem
// No resistance
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::CFDDEMProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;

// dem_mp
//  No resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;

// EPSD resistance
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  2,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
template class ParticleParticleContactForce<
  3,
  DEM::DEMMPProperties::PropertiesIndex,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::epsd_resistance>;
