#include <dem/dem_properties.h>
#include <dem/explicit_euler_integrator.h>

using namespace DEM;

template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
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

      // Calculate the acceleration of a particle
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::acc_x + d] =
            g[d] + particle_properties[PropertiesIndex::force_x + d] /
                     particle_properties[PropertiesIndex::mass];

          // Reinitializing force
          particle_properties[PropertiesIndex::force_x + d] = 0;
        }

      // Velocity integration:
      Tensor<1, dim> particle_velocity;

      for (unsigned int d = 0; d < dim; ++d)
        {
          particle_velocity[d] =
            particle_properties[PropertiesIndex::v_x + d] +
            dt * particle_properties[PropertiesIndex::acc_x + d];

          particle_properties[PropertiesIndex::v_x + d] = particle_velocity[d];
        }

      // Position integration
      auto particle_position = particle->get_location();
      particle_position      = particle_position + (particle_velocity * dt);
      particle->set_location(particle_position);

      // Angular velocity
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_properties[PropertiesIndex::M_x + d] /
                  particle_properties[PropertiesIndex::mom_inertia]);

          // Reinitializing torque
          particle_properties[PropertiesIndex::M_x + d] = 0;
        }
    }
}

template class ExplicitEulerIntegrator<2>;
template class ExplicitEulerIntegrator<3>;
