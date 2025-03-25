// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef multiphysics_integrator_h
#define multiphysics_integrator_h

// Deal.ii
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters_lagrangian.h>

#include <dem/dem_solver_parameters.h>

/**
 * @brief Implementation of an explicit euler scheme for the integration
 * of the particles' temperature.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param dt DEM time step.
 * @param torque Torque acting on particles.
 * @param force Force acting on particles.
 */
template <int dim, typename PropertiesIndex>
void
integrate_temperature(Particles::ParticleHandler<dim> &particle_handler,
                      const double                     dt,
                      std::vector<double>             &heat_transfer,
                      std::vector<double>             &heat_source);

#endif
