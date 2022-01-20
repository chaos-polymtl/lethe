#include <dem/dem_properties.h>
#include <dem/explicit_euler_integrator.h>

using namespace DEM;

// This function is empty for explicit Euler integrator
template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> & /*particle_handler*/,
  const Tensor<1, dim> & /*body_force*/,
  const double /*time_step*/,
  const std::vector<Tensor<1, dim>> & /*momentum*/,
  const std::vector<Tensor<1, dim>> & /*force*/,
  const std::vector<double> & /*MOI*/)
{}

template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, dim> &           g,
  const double                     dt,
  std::vector<Tensor<1, dim>> &    momentum,
  std::vector<Tensor<1, dim>> &    force,
  const std::vector<double> &      MOI)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties and location once
      // to improve efficiency
#if DEAL_II_VERSION_GTE(10, 0, 0)
      types::particle_index particle_id = particle->get_local_index();
#else
      types::particle_index particle_id = particle->get_id();
#endif
      auto            particle_properties = particle->get_properties();
      Tensor<1, dim> &particle_momentum   = momentum[particle_id];
      Tensor<1, dim> &particle_force      = force[particle_id];
      auto            particle_position   = particle->get_location();
      double mass_inverse = 1 / particle_properties[PropertiesIndex::mass];
      double MOI_inverse  = 1 / MOI[particle_id];


      for (int d = 0; d < dim; ++d)
        {
          acceleration[d] = g[d] + (particle_force[d]) * mass_inverse;

          // Velocity integration:
          particle_properties[PropertiesIndex::v_x + d] += dt * acceleration[d];

          // Reinitializing force
          particle_force[d] = 0;

          // Position integration
          particle_position[d] +=
            dt * particle_properties[PropertiesIndex::v_x + d];

          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_momentum[d] * MOI_inverse);

          // Reinitializing torque
          particle_momentum[d] = 0;
        }
      particle->set_location(particle_position);
    }
}

template class ExplicitEulerIntegrator<2>;
template class ExplicitEulerIntegrator<3>;
