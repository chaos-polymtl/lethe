#include <dem/dem_properties.h>
#include <dem/explicit_euler_integrator.h>

using namespace DEM;

// This function is empty for explicit Euler integrator
template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> & /*particle_handler*/,
  Tensor<1, dim> /*g*/,
  std::unordered_map<unsigned int, Tensor<1, dim>> & /*force*/,
  double /*dt*/)
{}

template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &                 particle_handler,
  Tensor<1, dim>                                    g,
  double                                            dt,
  std::unordered_map<unsigned int, Tensor<1, dim>> &momentum,
  std::unordered_map<unsigned int, Tensor<1, dim>> &force,
  std::unordered_map<unsigned int, double> &        MOI)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties and location once
      // to improve efficiency
      unsigned int    particle_id         = particle->get_id();
      auto            particle_MOI        = MOI[particle_id];
      auto            particle_properties = particle->get_properties();
      Tensor<1, dim> &particle_momentum   = momentum[particle_id];
      Tensor<1, dim> &particle_force      = force[particle_id];
      auto            particle_position   = particle->get_location();


      for (int d = 0; d < dim; ++d)
        {
          acceleration[d] = g[d] + particle_force[d] /
                                     particle_properties[PropertiesIndex::mass];

          // Velocity integration:
          particle_properties[PropertiesIndex::v_x + d] += dt * acceleration[d];

          // Reinitializing force
          particle_force[d] = 0;

          // Position integration
          particle_position[d] +=
            dt * particle_properties[PropertiesIndex::v_x + d];

          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_momentum[d] / particle_MOI);

          // Reinitializing torque
          particle_momentum[d] = 0;
        }
      particle->set_location(particle_position);
    }
}

template class ExplicitEulerIntegrator<2>;
template class ExplicitEulerIntegrator<3>;
