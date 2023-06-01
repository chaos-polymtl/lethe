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
  const std::vector<Tensor<1, 3>> & /*torque*/,
  const std::vector<Tensor<1, 3>> & /*force*/,
  const std::vector<double> & /*MOI*/)
{}

template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, 3> &             g,
  const double                     dt,
  std::vector<Tensor<1, 3>> &      torque,
  std::vector<Tensor<1, 3>> &      force,
  const std::vector<double> &      MOI)
{
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties and location once
      // to improve efficiency
#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
      types::particle_index particle_id = particle->get_id();
#else
      types::particle_index particle_id = particle->get_local_index();
#endif

      auto          particle_properties = particle->get_properties();
      Tensor<1, 3> &particle_torque     = torque[particle_id];
      Tensor<1, 3> &particle_force      = force[particle_id];
      Point<3>      particle_position;
      double mass_inverse = 1 / particle_properties[PropertiesIndex::mass];
      double MOI_inverse  = 1 / MOI[particle_id];


      if constexpr (dim == 3)
        particle_position = particle->get_location();

      if constexpr (dim == 2)
        particle_position = point_nd_to_3d(particle->get_location());

      for (int d = 0; d < 3; ++d)
        {
          acceleration[d] = g[d] + (particle_force[d]) * mass_inverse;

          // Velocity integration:
          particle_properties[PropertiesIndex::v_x + d] += dt * acceleration[d];

          // Position integration
          particle_position[d] +=
            dt * particle_properties[PropertiesIndex::v_x + d];

          particle_properties[PropertiesIndex::omega_x + d] +=
            dt * (particle_torque[d] * MOI_inverse);
        }

      // Reinitialize force
      particle_force = 0;

      // Reinitialize torque
      particle_torque = 0;

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

// Explicit Euler not implemented for disabling contacts
template <int dim>
void
ExplicitEulerIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> & /* particle_handler */,
  const Tensor<1, 3> & /* g */,
  const double /* dt */,
  std::vector<Tensor<1, 3>> & /* torque */,
  std::vector<Tensor<1, 3>> & /* force */,
  const std::vector<double> & /* MOI */,
  const parallel::distributed::Triangulation<dim> & /* triangulation */,
  DisableContacts<dim> & /* disable_contacts_object */)
{
  throw std::runtime_error(
    "Disabling particle contacts not supported with explicit Euler integrator, use Verlet integrator.");
}

template class ExplicitEulerIntegrator<2>;
template class ExplicitEulerIntegrator<3>;
