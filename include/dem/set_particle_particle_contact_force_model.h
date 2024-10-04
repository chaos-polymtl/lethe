// SPDX-FileCopyrightText: Copyright (c) 2022, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_set_particle_particle_contact_force_model_h
#define lethe_set_particle_particle_contact_force_model_h

#include <dem/dem_solver_parameters.h>
#include <dem/force_chains_visualization.h>
#include <dem/particle_particle_contact_force.h>


/**
 * @brief Return the particle-particle contact force model based on the
 * spring-dashpot model and the rolling resistance model (and the cohesive
 * force model if applicable).
 *
 * @tparam dim Space dimension.
 *
 * @param[in] dem_parameters DEM parameters.
 *
 * @return The pointer to the particle-particle contact force object.
 */
template <int dim>
std::shared_ptr<ParticleParticleContactForceBase<dim>>
set_particle_particle_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters);

/**
 *  @brief Set the rolling resistance model for the particle-particle contact
 *  force model.
 *
 * @tparam dim Space dimension.
 * @tparam particle_particle_contact_force_model Particle-particle contact force
 * model.
 *
 * @param[in] dem_parameters DEM parameters.
 * @param[out] particle_particle_contact_force_object Pointer to the particle-
 * particle contact force object.
 */
template <int dim,
          Parameters::Lagrangian::ParticleParticleContactForceModel
            particle_particle_contact_force_model>
void
set_rolling_resistance_model(
  const DEMSolverParameters<dim> &dem_parameters,
  std::shared_ptr<ParticleParticleContactForceBase<dim>>
    &particle_particle_contact_force_object);

/**
 * @brief Return the particle-particle contact force model based on the
 * spring-dashpot model and the rolling resistance model (and the cohesive
 * force model if applicable) for the force-chain post-processing.
 *
 * @param dem_parameters DEM parameters.
 *
 * @return The pointer to the particle-particle contact force object.
 */
template <int dim>
std::shared_ptr<ParticlesForceChainsBase<dim>>
set_force_chains_contact_force_model(
  const DEMSolverParameters<dim> &dem_parameters);

#endif
