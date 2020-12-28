#include <dem/dem_properties.h>
#include <dem/velocity_verlet_integrator.h>

using namespace DEM;

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate(
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
      auto           particle_properties = particle->get_properties();
      Tensor<1, dim> half_step_velocity;

      for (int d = 0; d < dim; ++d)
        {
          // Calculate the half step particle velocity
          half_step_velocity[d] =
            particle_properties[PropertiesIndex::v_x + d] +
            0.5 * dt * particle_properties[PropertiesIndex::acc_x + d];

          // Calculate the acceleration
          particle_properties[PropertiesIndex::acc_x + d] =
            g[d] + particle_properties[PropertiesIndex::force_x + d] /
                     particle_properties[PropertiesIndex::mass];

          // Reinitializing force
          particle_properties[PropertiesIndex::force_x + d] = 0;
        }

      // Calculate the particle position
      auto particle_position = particle->get_location();
      particle_position      = particle_position + (half_step_velocity * dt);
      particle->set_location(particle_position);


      // Update particle velocity
      for (int d = 0; d < dim; ++d)
        {
          // Calculate the particle full step velocity
          half_step_velocity[d] +=
            particle_properties[PropertiesIndex::acc_x + d] * 0.5 * dt;

          particle_properties[PropertiesIndex::v_x + d] = half_step_velocity[d];
        }

      // Store particle angular velocity and torque in Tensors
      Tensor<1, dim> particle_angular_velocity;
      Tensor<1, dim> particle_torque;
      for (int d = 0; d < dim; ++d)
        {
          particle_angular_velocity[d] =
            particle_properties[PropertiesIndex::omega_x + d];
          particle_torque[d] = particle_properties[PropertiesIndex::M_x + d];

          // Reinitializing torque
          particle_properties[PropertiesIndex::M_x + d] = 0;
        }

      // Updating angular velocity
      particle_angular_velocity =
        particle_angular_velocity +
        dt *
          (particle_torque / particle_properties[PropertiesIndex::mom_inertia]);

      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::omega_x + d] =
            particle_angular_velocity[d];
        }
    }
}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;
