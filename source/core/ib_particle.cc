#include <core/ib_particle.h>


template <int dim>
void
IBParticle<dim>::initialise_all()
{
  // initilise all the variables associated to an immersed boundary particle.
  mass               = 1;
  radius             = 1;
  local_alpha_torque = 1;
  local_alpha_force  = 1;
  particle_id        = 0;

  inertia[0][0] = 1;
  inertia[1][1] = 1;
  inertia[2][2] = 1;

  forces[0] = 0;
  forces[1] = 0;

  velocity[0] = 0;
  velocity[1] = 0;

  position[0] = 0;
  position[1] = 0;

  torques[0] = 0;
  torques[1] = 0;
  torques[2] = 0;

  angular_position[0] = 0;
  angular_position[1] = 0;
  angular_position[2] = 0;

  omega[0] = 0;
  omega[1] = 0;
  omega[2] = 0;

  if (dim == 3)
    {
      forces[2]   = 0;
      velocity[2] = 0;
      position[2] = 0;
    }

  last_forces           = forces;
  last_torques          = torques;
  last_position         = position;
  last_velocity         = velocity;
  velocity_iter         = velocity;
  last_omega            = omega;
  omega_iter            = omega;
  omega_impulsion       = 0;
  omega_impulsion_iter  = 0;
  last_angular_position = angular_position;
  impulsion             = 0;
  impulsion_iter        = 0;


  youngs_modulus=10000000;
  restitution_coefficient=1;
  friction_coefficient=1;
  poisson_ratio=0.3;
  rolling_friction_coefficient=0.1;
}

template <int dim>
void
IBParticle<dim>::initialise_last()
{
  // initilise all the variables associated to an immersed boundary particle
  last_forces           = forces;
  last_position         = position;
  last_velocity         = velocity;
  velocity_iter         = velocity;
  impulsion_iter        = impulsion;
  last_omega            = omega;
  omega_iter            = omega;
  omega_impulsion_iter  = omega_impulsion;
  last_angular_position = angular_position;
}

template <int dim>
std::vector<std::pair<std::string, int>>
IBParticle<dim>::get_properties_name()
{
  std::vector<std::pair<std::string, int>> properties(
    PropertiesIndex::n_properties);
  properties[PropertiesIndex::id] = std::make_pair("ID", 1);
  properties[PropertiesIndex::dp] = std::make_pair("Diameter", 1);
  properties[PropertiesIndex::vx] = std::make_pair("Velocity", 3);
  properties[PropertiesIndex::vy] = std::make_pair("Velocity", 1);
  properties[PropertiesIndex::vz] = std::make_pair("Velocity", 1);
  properties[PropertiesIndex::fx] = std::make_pair("Force", 3);
  properties[PropertiesIndex::fy] = std::make_pair("Force", 1);
  properties[PropertiesIndex::fz] = std::make_pair("Force", 1);
  properties[PropertiesIndex::ox] = std::make_pair("Omega", 3);
  properties[PropertiesIndex::oy] = std::make_pair("Omega", 1);
  properties[PropertiesIndex::oz] = std::make_pair("Omega", 1);
  properties[PropertiesIndex::tx] = std::make_pair("Torques", 3);
  properties[PropertiesIndex::ty] = std::make_pair("Torques", 1);
  properties[PropertiesIndex::tz] = std::make_pair("Torques", 1);

  return properties;
}
template <int dim>
std::vector<double>
IBParticle<dim>::get_properties()
{
  std::vector<double> properties(get_number_properties());
  properties[0]  = particle_id;
  properties[1]  = radius * 2;
  properties[2]  = velocity[0];
  properties[3]  = velocity[1];
  properties[5]  = forces[0];
  properties[6]  = forces[1];
  properties[8]  = omega[0];
  properties[9]  = omega[1];
  properties[10] = omega[2];
  properties[11] = torques[0];
  properties[12] = torques[1];
  properties[13] = torques[2];
  if (dim == 2)
    {
      properties[4] = 0;
      properties[7] = 0;
    }
  if (dim == 3)
    {
      properties[4] = velocity[2];
      properties[7] = forces[2];
    }

  return properties;
}

template <int dim>
unsigned int
IBParticle<dim>::get_number_properties()
{
  return PropertiesIndex::n_properties;
}

template class IBParticle<2>;
template class IBParticle<3>;
