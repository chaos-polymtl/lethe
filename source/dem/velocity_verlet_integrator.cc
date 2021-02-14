#include <dem/dem_properties.h>
#include <dem/velocity_verlet_integrator.h>

using namespace DEM;

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate_pre_force(
  Particles::ParticleHandler<dim> &particle_handler,
  Tensor<1, dim> /*g*/,
  double dt)
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
  Particles::ParticleHandler<dim> &        particle_handler,
  Tensor<1, dim>                           g,
  double                                   dt,
  std::unordered_map<int, Tensor<1, dim>> &momentum,
  std::unordered_map<int, Tensor<1, dim>> &force,
  std::unordered_map<int, double> &        MOI,
  std::unordered_map<int, Tensor<1, dim>> &acceleration)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      unsigned int    particle_id           = particle->get_id();
      auto            particle_MOI          = MOI[particle_id];
      auto            particle_properties   = particle->get_properties();
      Tensor<1, dim> &particle_momentum     = momentum[particle_id];
      Tensor<1, dim> &particle_force        = force[particle_id];
      Tensor<1, dim> &particle_acceleration = acceleration[particle_id];

      for (int d = 0; d < dim; ++d)
        {
          // Calculate the acceleration
          particle_acceleration[d] =
            g[d] +
            particle_force[d] / particle_properties[PropertiesIndex::mass];

          // Reinitializing force
          particle_force[d] = 0;

          // Calculate the particle full step velocity
          particle_properties[PropertiesIndex::v_x + d] +=
            particle_acceleration[d] * 0.5 * dt;

          // Updating angular velocity
          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_momentum[d] / particle_MOI);

          // Reinitializing torque
          particle_momentum[d] = 0;

          // Calculate the half step particle velocity for the next integration
          // (pre-force) step. I moved this part from the beginning of pre-force
          // integration.
          particle_properties[PropertiesIndex::v_x + d] +=
            0.5 * dt * particle_acceleration[d];
        }
    }
}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;
