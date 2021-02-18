#include <dem/dem_properties.h>
#include <dem/velocity_verlet_integrator.h>

using namespace DEM;

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> &                          particle_handler,
  Tensor<1, dim> &                                           g,
  std::unordered_map<types::particle_index, Tensor<1, dim>> &force,
  double                                                     dt)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      auto            particle_properties = particle->get_properties();
      auto            particle_position   = particle->get_location();
      unsigned int    particle_id         = particle->get_id();
      Tensor<1, dim>  particle_acceleration;
      Tensor<1, dim> &particle_force = force[particle_id];

      for (int d = 0; d < dim; ++d)
        {
          // Update acceleration
          particle_acceleration[d] =
            g[d] +
            particle_force[d] / particle_properties[PropertiesIndex::mass];

          // Half-step velocity
          particle_properties[PropertiesIndex::v_x + d] +=
            0.5 * particle_acceleration[d] * dt;

          // Update particle position using half-step velocity
          particle_position[d] +=
            particle_properties[PropertiesIndex::v_x + d] * dt;
        }
      particle->set_location(particle_position);
    }
}

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &                          particle_handler,
  Tensor<1, dim> &                                           g,
  double                                                     dt,
  std::unordered_map<types::particle_index, Tensor<1, dim>> &momentum,
  std::unordered_map<types::particle_index, Tensor<1, dim>> &force,
  std::unordered_map<types::particle_index, double> &        MOI)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      unsigned int    particle_id         = particle->get_id();
      auto            particle_properties = particle->get_properties();
      Tensor<1, dim> &particle_momentum   = momentum[particle_id];
      Tensor<1, dim> &particle_force      = force[particle_id];
      Tensor<1, dim>  particle_acceleration;
      auto            particle_position = particle->get_location();
      double mass_inverse = 1 / particle_properties[PropertiesIndex::mass];
      double MOI_inverse  = 1 / MOI[particle_id];

      for (int d = 0; d < dim; ++d)
        {
          particle_acceleration[d] = g[d] + particle_force[d] * mass_inverse;

          // Particle velocity integration
          particle_properties[PropertiesIndex::v_x + d] +=
            dt * particle_acceleration[d];

          // Particle location integration
          particle_position[d] +=
            particle_properties[PropertiesIndex::v_x + d] * dt;

          // Updating angular velocity
          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_momentum[d] * MOI_inverse);

          // Reinitializing force
          particle_force[d] = 0;

          // Reinitializing torque
          particle_momentum[d] = 0;
        }
      particle->set_location(particle_position);
    }
}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;
