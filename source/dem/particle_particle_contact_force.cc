// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_lagrangian.h>

#include <dem/particle_particle_contact_force.h>

using namespace DEM;
using namespace Parameters::Lagrangian;

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
{
  set_effective_properties(dem_parameters);
}

template <int                               dim,
          ParticleParticleContactForceModel contact_model,
          RollingResistanceMethod           rolling_friction_model>
void
ParticleParticleContactForce<dim, contact_model, rolling_friction_model>::
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
  // Calculating the contact forces the local-local adjacent particles.
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::local_particle_particle>(
        adjacent_particles_list, torque, force, dt);
    }

  // Calculating the contact forces the local-ghost adjacent particles.
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::ghost_particle_particle>(
        adjacent_particles_list, torque, force, dt);
    }

  // Calculating the contact forces the local-local periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::local_periodic_particle_particle>(
        periodic_adjacent_particles_list, torque, force, dt);
    }

  // Calculating the contact forces the local-ghost periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_ghost_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_periodic_particle_particle>(
        periodic_adjacent_particles_list, torque, force, dt);
    }

  // Calculating the contact forces the ghost-local periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       ghost_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_local_periodic_particle_particle>(
        periodic_adjacent_particles_list, torque, force, dt);
    }
}

// No resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::no_resistance>;

// Constant resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::constant_resistance>;

// Viscous resistance
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::DMT,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_JKR,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_force,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::hertz_mindlin_limit_overlap,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  2,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
template class ParticleParticleContactForce<
  3,
  ParticleParticleContactForceModel::linear,
  RollingResistanceMethod::viscous_resistance>;
