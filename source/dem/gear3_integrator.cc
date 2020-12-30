#include <dem/dem_properties.h>
#include <dem/gear3_integrator.h>

using namespace DEM;

template <int dim>
void
Gear3Integrator<dim>::integrate(
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
      auto particle_position   = particle->get_location();

      // Updating particle location, velocity, acceleration and derivative of
      // acceleration:
      for (int d = 0; d < dim; ++d)
        {
          // Predictor
          predicted_location[d] =
            particle_position[d] +
            (particle_properties[PropertiesIndex::v_x + d] * dt) +
            (particle_properties[PropertiesIndex::acc_x + d] * dt * dt * 0.5) +
            (particle_properties[PropertiesIndex::acc_derivative_x + d] * dt *
             dt * dt * 0.1667);
          predicted_velocity[d] =
            particle_properties[PropertiesIndex::v_x + d] +
            (particle_properties[PropertiesIndex::acc_x + d] * dt) +
            (particle_properties[PropertiesIndex::acc_derivative_x + d] * dt *
             dt * 0.5);
          predicted_acceleration[d] =
            particle_properties[PropertiesIndex::acc_x + d] +
            (particle_properties[PropertiesIndex::acc_derivative_x + d] * dt);

          // Finding corrected acceleration
          corrected_accereration[d] =
            g[d] + particle_properties[PropertiesIndex::force_x + d] /
                     particle_properties[PropertiesIndex::mass];

          // Reinitializing force
          particle_properties[PropertiesIndex::force_x + d] = 0;

          // Calculation of acceleration deviation
          acceleration_deviation[d] =
            corrected_accereration[d] - predicted_acceleration[d];

          // Corrector
          corrected_location[d] =
            predicted_location[d] +
            acceleration_deviation[d] * (0.0833 * dt * dt);
          particle_properties[PropertiesIndex::v_x + d] =
            predicted_velocity[d] + acceleration_deviation[d] * (0.4167 * dt);
          particle_properties[PropertiesIndex::acc_x + d] =
            predicted_acceleration[d] + acceleration_deviation[d];
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
}

template class Gear3Integrator<2>;
template class Gear3Integrator<3>;
