// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef integrate_temperature_h
#define integrate_temperature_h

#include <deal.II/particles/particle_handler.h>


using namespace DEMMP;

template <int dim, typename PropertiesIndex>
void
integrate_temperature(Particles::ParticleHandler<dim> &particle_handler,
                      const double                     dt,
                      std::vector<double>             &heat_transfer,
                      const double                     heat_source)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      types::particle_index particle_id = particle->get_local_index();

      auto    particle_properties    = particle->get_properties();
      double &particle_heat_transfer = heat_transfer[particle_id];

      double mass_inverse = 1 / particle_properties[PropertiesIndex::mass];
      double heat_capacity_inverse =
        1 / particle_properties[PropertiesIndex::heat_capacity];

      // Integration
      particle_properties[PropertiesIndex::T] +=
        dt * (particle_heat_transfer + heat_source) * mass_inverse *
        heat_capacity_inverse;

      // Reinitialize heat_transfer
      particle_heat_transfer = 0;
    }
}


#endif
