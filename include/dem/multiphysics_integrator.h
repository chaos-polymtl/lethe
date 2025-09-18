// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef multiphysics_integrator_h
#define multiphysics_integrator_h

// Deal.ii
#include <dem/dem_solver_parameters.h>

#include <deal.II/particles/particle_handler.h>

/**
 * @brief Implementation of an explicit euler time-stepping scheme for the integration
 * of the particles' temperature.
 *
 * @tparam dim Integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 * @param particle_handler Storage of particles and their accessor functions.
 * @param[in] dt DEM time step.
 * @param[in,out] heat_transfer_rate Particle-particle heat transfer rate.
 * @param[in] heat_source Additional heat source term.
 */
template <int dim, typename PropertiesIndex>
void
integrate_temperature(Particles::ParticleHandler<dim> &particle_handler,
                      const double                     dt,
                      std::vector<double>             &heat_transfer_rate,
                      const std::vector<double>       &heat_source);

#endif
