// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef integrate_temperature_h
#define integrate_temperature_h

// Deal.ii
#include <deal.II/particles/particle_handler.h>

// Lethe
#include <core/dem_properties.h>
#include <core/parameters_lagrangian.h>

#include <dem/dem_solver_parameters.h>


template <int dim, typename PropertiesIndex>
void
integrate_temperature(Particles::ParticleHandler<dim> &particle_handler,
                      const DEMSolverParameters<dim>  &dem_parameters,
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
      unsigned int type   = particle_properties[PropertiesIndex::type];
      double       heat_capacity =
        dem_parameters.lagrangian_physical_properties.heat_capacity_particle.at(
          type);

      // Integration
      particle_properties[PropertiesIndex::T] +=
        dt * (particle_heat_transfer + heat_source) * mass_inverse * 1 /
        heat_capacity;

      // Reinitialize heat_transfer
      particle_heat_transfer = 0;
    }
}

#endif
