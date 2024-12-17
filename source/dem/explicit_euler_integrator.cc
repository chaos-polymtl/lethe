// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/explicit_euler_integrator.h>

using namespace DEM;

// This function is empty for explicit Euler integrator
template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> & /*particle_handler*/,
  const Tensor<1, 3> & /*body_force*/,
  const double /*time_step*/,
  const std::vector<double> & /*MOI*/)
{}

template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, 3>              &g,
  const double                     dt,
  const std::vector<double>       &MOI)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties and location once
      // to improve efficiency
      types::particle_index particle_id = particle->get_local_index();

      auto particle_properties = particle->get_properties();

      Point<3> particle_position;
      double   mass_inverse = 1 / particle_properties[PropertiesIndex::mass];
      double   MOI_inverse  = 1 / MOI[particle_id];


      if constexpr (dim == 3)
        particle_position = particle->get_location();

      if constexpr (dim == 2)
        particle_position = point_nd_to_3d(particle->get_location());

      Tensor<1, 3> acceleration;
      for (int d = 0; d < 3; ++d)
        {
          // Velocity integration:
          particle_properties[PropertiesIndex::v_x + d] +=
            dt * (g[d] + (particle_properties[PropertiesIndex::force_x + d]) *
                           mass_inverse);

          // Position integration
          particle_position[d] +=
            dt * particle_properties[PropertiesIndex::v_x + d];

          particle_properties[PropertiesIndex::omega_x + d] +=
            dt *
            (particle_properties[PropertiesIndex::torque_x + d] * MOI_inverse);

          // Reset torque and force to zero for the next iteration
          particle_properties[PropertiesIndex::torque_x + d] = 0;
          particle_properties[PropertiesIndex::force_x + d]  = 0;
        }

      if constexpr (dim == 3)
        particle->set_location(particle_position);

      if constexpr (dim == 2)
        {
          Point<2> position_2d;
          position_2d[0] = particle_position[0];
          position_2d[1] = particle_position[1];
          particle->set_location(position_2d);
        }
    }
}

// Explicit Euler not implemented for adaptive sparse contacts
template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, 3>              &g,
  const double                     dt,
  const std::vector<double>       &MOI,
  const parallel::distributed::Triangulation<dim> & /* triangulation */,
  AdaptiveSparseContacts<dim> & /* sparse_contacts_object */)
{
  auto *action_manager = DEMActionManager::get_action_manager();

  bool use_default_function =
    !action_manager->check_sparse_contacts_enabled() ||
    action_manager->check_mobility_status_reset();

  if (use_default_function)
    {
      integrate(particle_handler, g, dt, MOI);
      return;
    }

  throw std::runtime_error(
    "Adaptive sparse contacts are not supported with explicit Euler integrator, use Velocity Verlet integrator.");
}

template class ExplicitEulerIntegrator<2>;
template class ExplicitEulerIntegrator<3>;
