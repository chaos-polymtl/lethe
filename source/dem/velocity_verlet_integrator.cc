#include <core/dem_properties.h>
#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/velocity_verlet_integrator.h>

using namespace DEM;

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate_half_step_location(
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, 3> &             g,
  const double                     dt,
  const std::vector<Tensor<1, 3>> &torque,
  const std::vector<Tensor<1, 3>> &force,
  const std::vector<double> &      MOI)
{
  Tensor<1, dim> particle_acceleration;

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      auto                  particle_properties = particle->get_properties();
      auto                  particle_position   = particle->get_location();
      types::particle_index particle_id         = particle->get_local_index();

      for (int d = 0; d < dim; ++d)
        {
          // Update acceleration
          particle_acceleration[d] =
            g[d] + (force[particle_id][d]) /
                     particle_properties[PropertiesIndex::mass];

          // Half-step velocity
          particle_properties[PropertiesIndex::v_x + d] +=
            0.5 * particle_acceleration[d] * dt;

          // Half-step angular velocity
          particle_properties[PropertiesIndex::omega_x + d] +=
            0.5 * (torque[particle_id][d] / MOI[particle_id]) * dt;

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
  Particles::ParticleHandler<dim> &particle_handler,
  const Tensor<1, 3> &             g,
  const double                     dt,
  std::vector<Tensor<1, 3>> &      torque,
  std::vector<Tensor<1, 3>> &      force,
  const std::vector<double> &      MOI)
{
  Point<3>           particle_position;
  const Tensor<1, 3> dt_g = g * dt;

  for (auto &particle : particle_handler)
    {
      // Get the total array view to the particle properties once to improve
      // efficiency
      types::particle_index particle_id = particle.get_local_index();

      auto          particle_properties = particle.get_properties();
      Tensor<1, 3> &particle_torque     = torque[particle_id];
      Tensor<1, 3> &particle_force      = force[particle_id];

      double dt_mass_inverse = dt / particle_properties[PropertiesIndex::mass];
      double dt_MOI_inverse  = dt / MOI[particle_id];

      particle_position = [&] {
        if constexpr (dim == 3)
          {
            return particle.get_location();
          }
        else
          {
            return (point_nd_to_3d(particle.get_location()));
          }
      }();

      // Loop is manually unrolled for performance reason.
      // for (int d = 0; d < 3; ++d)
      {
        // Particle velocity integration
        particle_properties[PropertiesIndex::v_x] +=
          dt_g[0] + particle_force[0] * dt_mass_inverse;
        particle_properties[PropertiesIndex::v_y] +=
          dt_g[1] + particle_force[1] * dt_mass_inverse;
        particle_properties[PropertiesIndex::v_z] +=
          dt_g[2] + particle_force[2] * dt_mass_inverse;


        // Particle location integration
        particle_position[0] += particle_properties[PropertiesIndex::v_x] * dt;
        particle_position[1] += particle_properties[PropertiesIndex::v_y] * dt;
        particle_position[2] += particle_properties[PropertiesIndex::v_z] * dt;

        // Updating angular velocity
        particle_properties[PropertiesIndex::omega_x] +=
          particle_torque[0] * dt_MOI_inverse;
        particle_properties[PropertiesIndex::omega_y] +=
          particle_torque[1] * dt_MOI_inverse;
        particle_properties[PropertiesIndex::omega_z] +=
          particle_torque[2] * dt_MOI_inverse;
      }

      // Reinitialize force and torque of particle
      particle_force  = 0;
      particle_torque = 0;

      if constexpr (dim == 3)
        particle.set_location(particle_position);

      if constexpr (dim == 2)
        {
          Point<2> position_2d;
          position_2d[0] = particle_position[0];
          position_2d[1] = particle_position[1];
          particle.set_location(position_2d);
        }
    }
}

template <int dim>
void
VelocityVerletIntegrator<dim>::integrate(
  Particles::ParticleHandler<dim> &                particle_handler,
  const Tensor<1, 3> &                             g,
  const double                                     dt,
  std::vector<Tensor<1, 3>> &                      torque,
  std::vector<Tensor<1, 3>> &                      force,
  const std::vector<double> &                      MOI,
  const parallel::distributed::Triangulation<dim> &triangulation,
  std::unordered_map<types::global_cell_index, unsigned int>
    &cell_mobility_status_map)
{
  Point<3>           particle_position;
  const Tensor<1, 3> dt_g = g * dt;

  for (auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Check if the cell is mobile, if not, integration is skipped
          unsigned int mobile_status  = 2; // Value of the mobile status
          bool         cell_is_mobile = false;

          // Get the iterator to the cell in the map
          auto it = cell_mobility_status_map.find(cell->active_cell_index());

          // Check if the cell exist in the map (meaning that the cell is mobile
          // or active, inactive cells are not stored) and verify if the status
          // is mobile : true = iterator exist && value is mobile
          cell_is_mobile =
            it != cell_mobility_status_map.end() && it->second == mobile_status;

          // We loop over the particles, even if cell is not mobile, to reset
          // the force and torques value of the particle with the particle id
          auto particles_in_cell = particle_handler.particles_in_cell(cell);
          for (auto &particle : particles_in_cell)
            {
              types::particle_index particle_id = particle.get_local_index();
              if (cell_is_mobile)
                {
                  // Get the total array view to the particle properties once to
                  // improve efficiency
                  auto particle_properties = particle.get_properties();

                  Tensor<1, 3> &particle_torque = torque[particle_id];
                  Tensor<1, 3> &particle_force  = force[particle_id];

                  double dt_mass_inverse =
                    dt / particle_properties[PropertiesIndex::mass];
                  double dt_MOI_inverse = dt / MOI[particle_id];

                  particle_position = [&] {
                    if constexpr (dim == 3)
                      {
                        return particle.get_location();
                      }
                    else
                      {
                        return (point_nd_to_3d(particle.get_location()));
                      }
                  }();

                  // Loop is manually unrolled for performance reason.
                  // for (int d = 0; d < 3; ++d)
                  {
                    // Particle velocity integration
                    particle_properties[PropertiesIndex::v_x] +=
                      dt_g[0] + particle_force[0] * dt_mass_inverse;
                    particle_properties[PropertiesIndex::v_y] +=
                      dt_g[1] + particle_force[1] * dt_mass_inverse;
                    particle_properties[PropertiesIndex::v_z] +=
                      dt_g[2] + particle_force[2] * dt_mass_inverse;


                    // Particle location integration
                    particle_position[0] +=
                      particle_properties[PropertiesIndex::v_x] * dt;
                    particle_position[1] +=
                      particle_properties[PropertiesIndex::v_y] * dt;
                    particle_position[2] +=
                      particle_properties[PropertiesIndex::v_z] * dt;

                    // Updating angular velocity
                    particle_properties[PropertiesIndex::omega_x] +=
                      particle_torque[0] * dt_MOI_inverse;
                    particle_properties[PropertiesIndex::omega_y] +=
                      particle_torque[1] * dt_MOI_inverse;
                    particle_properties[PropertiesIndex::omega_z] +=
                      particle_torque[2] * dt_MOI_inverse;
                  }

                  // Reinitialize force and torque of particle
                  particle_force  = 0;
                  particle_torque = 0;

                  if constexpr (dim == 3)
                    particle.set_location(particle_position);

                  if constexpr (dim == 2)
                    {
                      Point<2> position_2d;
                      position_2d[0] = particle_position[0];
                      position_2d[1] = particle_position[1];
                      particle.set_location(position_2d);
                    }
                }
              else
                {
                  // Reset forces
                  torque[particle_id] = 0.0;
                  force[particle_id]  = 0.0;
                }
            }
        }
    }
}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;
