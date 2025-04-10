// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/multiphysics_integrator.h>

template <int dim, typename PropertiesIndex>
void
integrate_temperature(Particles::ParticleHandler<dim> &particle_handler,
                      const double                     dt,
                      std::vector<double>             &heat_transfer_rate,
                      const std::vector<double>       &heat_source)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      const types::particle_index particle_id = particle->get_local_index();
      auto         particle_properties        = particle->get_properties();
      double      &particle_heat_transfer_rate     = heat_transfer_rate[particle_id];
      double       particle_heat_source       = heat_source[particle_id];
      const double mass_inverse =
        1 / particle_properties[PropertiesIndex::mass];
      const double specific_heat_inverse =
        1 / particle_properties[PropertiesIndex::specific_heat];

      // Integration
      particle_properties[PropertiesIndex::T] +=
        dt * (particle_heat_transfer_rate + particle_heat_source) * mass_inverse *
        specific_heat_inverse;

      // Reinitialize heat_transfer_rate
      particle_heat_transfer_rate = 0;
    }
}

template void
integrate_temperature<2, DEM::DEMMPProperties::PropertiesIndex>(
  Particles::ParticleHandler<2> &particle_handler,
  const double                   dt,
  std::vector<double>           &heat_transfer_rate,
  const std::vector<double>     &heat_source);

template void
integrate_temperature<3, DEM::DEMMPProperties::PropertiesIndex>(
  Particles::ParticleHandler<3> &particle_handler,
  const double                   dt,
  std::vector<double>           &heat_transfer_rate,
  const std::vector<double>     &heat_source);
