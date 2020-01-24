#include <dem/particle_properties_index.h>
#include <dem/velocity_verlet_integrator.h>


using namespace DEM;

template <int dim, int spacedim>
void
VelocityVerletIntegrator<dim, spacedim>::integrate(
  Particles::ParticleHandler<dim, spacedim> &particle_handler,
  Tensor<1, dim>                             g,
  double                                     dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      auto particle_properties = particle->get_properties();

      // Calculate the acceleration of a particle
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::acc_x + d] =
            g[d] + particle_properties[PropertiesIndex::force_x + d] /
                     particle_properties[PropertiesIndex::mass];
        }

      // Store particle velocity and acceleration in Tensors
      Tensor<1, dim> particle_velocity;
      Tensor<1, dim> particle_acceleration;
      for (int d = 0; d < dim; ++d)
        {
          particle_velocity[d] = particle_properties[PropertiesIndex::v_x + d];
          particle_acceleration[d] =
            particle_properties[PropertiesIndex::acc_x + d];
        }

      // Calculate the half step particle velocity
      Tensor<1, dim> vStar =
        particle_velocity + (particle_acceleration * 0.5 * dt);

      // Calculate the particle position
      auto particle_position = particle->get_location();
      particle_position      = particle_position + (particle_velocity * dt);
      particle->set_location(particle_position);

      // Calculate the particle full step velocity
      vStar = vStar + particle_acceleration * 0.5 * dt;

      // Update particle velocity
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::v_x + d] = vStar[d];
        }

      // Angular velocity:
      /*
      particle->get_properties()[16] = particle->get_properties()[16] +
        (particle->get_properties()[21]) / (particle->get_properties()[20]);
      particle->get_properties()[17] =particle->get_properties()[17] +
        (particle->get_properties()[22]) / (particle->get_properties()[20]);
      particle->get_properties()[18] = particle->get_properties()[18] +
        (particle->get_properties()[23]) / (particle->get_properties()[20]);
        */
    }
}

template class VelocityVerletIntegrator<3, 3>;
