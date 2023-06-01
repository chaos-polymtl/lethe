#include <core/dem_properties.h>

#include <dem/gear3_integrator.h>

using namespace DEM;

template <int dim>
void
Gear3Integrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> & /*particle_handler*/,
  const Tensor<1, 3> & /*body_force*/,
  const double /*time_step*/,
  const std::vector<Tensor<1, 3>> & /*torque*/,
  const std::vector<Tensor<1, 3>> & /*force*/,
  const std::vector<double> & /*MOI*/)
{}

template <int dim>
void
Gear3Integrator<dim>::integrate(
  Particles::ParticleHandler<dim> & /*particle_handler*/,
  const Tensor<1, 3> & /*g*/,
  const double /*dt*/,
  std::vector<Tensor<1, 3>> & /*torque*/,
  std::vector<Tensor<1, 3>> & /*force*/,
  const std::vector<double> & /*MOI*/)
{
  /*
for (auto particle = particle_handler.begin();
     particle != particle_handler.end();
     ++particle)
  {
    // Get the total array view to the particle properties once to improve
    // efficiency
    auto particle_properties = particle->get_properties();

    // Updating particle location, velocity, acceleration and derivative of
    // acceleration:
    for (int d = 0; d < dim; ++d)
      {
        // Finding corrected acceleration
        corrected_accereration[d] =
          g[d] + particle_properties[PropertiesIndex::force_x + d] /
                   particle_properties[PropertiesIndex::mass];

        // Reinitializing force
        particle_properties[PropertiesIndex::force_x + d] = 0;

        // Calculation of acceleration deviation
        acceleration_deviation[d] =
          corrected_accereration[d] -
          particle_properties[PropertiesIndex::acc_x + d];

        // Corrector
        corrected_location[d] =
          predicted_location[d] +
          acceleration_deviation[d] * (0.0833 * dt * dt);
        particle_properties[PropertiesIndex::v_x + d] =
          particle_properties[PropertiesIndex::v_x + d] +
          acceleration_deviation[d] * (0.4167 * dt);
        particle_properties[PropertiesIndex::acc_x + d] =
          particle_properties[PropertiesIndex::v_x + d] +
          acceleration_deviation[d];
        particle_properties[PropertiesIndex::acc_derivative_x + d] +=
          acceleration_deviation[d] / dt;

        // Angular velocity
        particle_properties[PropertiesIndex::omega_x + d] +=
          dt * (particle_properties[PropertiesIndex::M_x + d] /
                particle_properties[PropertiesIndex::mom_inertia]);

        // Reinitializing torque
        particle_properties[PropertiesIndex::M_x + d] = 0;
      }
    particle->set_location(corrected_location);
  }
  */
}

// Gear 3 not implemented for disabling contacts
template <int dim>
void
Gear3Integrator<dim>::integrate(
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
    "Disabling particle contacts not supported with explicit Gear 3 integrator, use Verlet integrator.");
}

template class Gear3Integrator<2>;
template class Gear3Integrator<3>;
