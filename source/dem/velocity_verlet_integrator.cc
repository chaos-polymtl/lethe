#include <dem/dem_properties.h>
#include <dem/velocity_verlet_integrator.h>

using namespace DEM;

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate_pre_force(
  Particles::ParticleHandler<dim> &particle_handler,
  double                           dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      auto particle_properties = particle->get_properties();
      auto particle_position   = particle->get_location();

      for (int d = 0; d < dim; ++d)
        {
          // Calculate the half step particle velocity
          particle_properties[PropertiesIndex::v_x + d] +=
            0.5 * dt * particle_properties[PropertiesIndex::acc_x + d];

          // Update particle position
          particle_position[d] +=
            (particle_properties[PropertiesIndex::v_x + d] * dt);
        }
      particle->set_location(particle_position);
    }
}

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate_post_force(
  Particles::ParticleHandler<dim> &particle_handler,
  Tensor<1, dim>                   g,
  double                           dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      auto particle_properties = particle->get_properties();

      for (int d = 0; d < dim; ++d)
        {
          // Calculate the acceleration
          particle_properties[PropertiesIndex::acc_x + d] =
            g[d] + particle_properties[PropertiesIndex::force_x + d] /
                     particle_properties[PropertiesIndex::mass];

          // Reinitializing force
          particle_properties[PropertiesIndex::force_x + d] = 0;

          // Calculate the particle full step velocity
          particle_properties[PropertiesIndex::v_x + d] +=
            particle_properties[PropertiesIndex::acc_x + d] * 0.5 * dt;

          // Updating angular velocity
          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_properties[PropertiesIndex::M_x + d] /
                  particle_properties[PropertiesIndex::mom_inertia]);

          // Reinitializing torque
          particle_properties[PropertiesIndex::M_x + d] = 0;
        }
    }
}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;
