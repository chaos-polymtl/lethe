// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_contact_force_templates_h
#define lethe_particle_particle_contact_force_templates_h

// This file contains the definitions of the member functions of
// ParticleParticleContactForce that are not defined in the class body.
// They are in a header rather than in a source file because the explicit
// instantiations of the class are distributed over two source files
// (particle_particle_contact_force_2d.cc and
// particle_particle_contact_force_3d.cc) so that they are compiled in parallel.
// Only these source files should include this header; everything else should
// include particle_particle_contact_force.h.

#include <core/parameters_lagrangian.h>

#include <dem/particle_particle_contact_force.h>

template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
ParticleParticleContactForce<dim,
                             PropertiesIndex,
                             contact_model,
                             rolling_friction_model>::
  ParticleParticleContactForce(const DEMSolverParameters<dim> &dem_parameters)
  : dmt_cut_off_threshold(dem_parameters.model_parameters.dmt_cut_off_threshold)
  , f_coefficient_epsd(dem_parameters.model_parameters.f_coefficient_epsd)
{
  set_effective_properties(dem_parameters);
  if constexpr (DEM::has_thermal_properties<PropertiesIndex>)
    {
      set_multiphysic_properties(dem_parameters);
    }
}

template <
  int dim,
  typename PropertiesIndex,
  Parameters::Lagrangian::ParticleParticleContactForceModel contact_model,
  Parameters::Lagrangian::RollingResistanceMethod rolling_friction_model>
void
ParticleParticleContactForce<dim,
                             PropertiesIndex,
                             contact_model,
                             rolling_friction_model>::
  calculate_particle_particle_contact(
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &local_adjacent_particles,
    typename DEM::dem_data_structures<dim>::adjacent_particle_pairs
      &ghost_adjacent_particles,
    typename DEM::dem_data_structures<dim>::periodic_adjacent_particle_pairs
      &local_local_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::periodic_adjacent_particle_pairs
      &local_ghost_periodic_adjacent_particles,
    typename DEM::dem_data_structures<dim>::periodic_adjacent_particle_pairs
                &ghost_local_periodic_adjacent_particles,
    const double dt,
    ParticleInteractionOutcomes<PropertiesIndex> &contact_outcome)
{
  // Calculating the contact forces and heat transfer rates for local-local
  // adjacent particles.
  for (auto &&adjacent_particles_list :
       local_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::local_particle_particle>(
        adjacent_particles_list, dt, contact_outcome);
    }

  // Calculating the contact forces and heat transfer rates for local-ghost
  // adjacent particles.
  for (auto &&adjacent_particles_list :
       ghost_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<ContactType::ghost_particle_particle>(
        adjacent_particles_list, dt, contact_outcome);
    }

  // Calculating the contact forces and heat transfer rates for local-local
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::local_periodic_particle_particle>(
        periodic_adjacent_particles_list, dt, contact_outcome);
    }

  // Calculating the contact forces and heat transfer rates for local-ghost
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       local_ghost_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_periodic_particle_particle>(
        periodic_adjacent_particles_list, dt, contact_outcome);
    }

  // Calculating the contact forces and heat transfer rates for ghost-local
  // periodic adjacent particles.
  for (auto &&periodic_adjacent_particles_list :
       ghost_local_periodic_adjacent_particles | boost::adaptors::map_values)
    {
      execute_contact_calculation<
        ContactType::ghost_local_periodic_particle_particle>(
        periodic_adjacent_particles_list, dt, contact_outcome);
    }
}

#endif
