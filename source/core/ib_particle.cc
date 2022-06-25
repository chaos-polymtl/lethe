#include <core/ib_particle.h>
#include <core/shape.h>


template <int dim>
void
IBParticle<dim>::initialise_all()
{
  // initialise all the variables associated to an immersed boundary particle.
  radius      = 1;
  particle_id = 0;

  inertia[0][0] = 1;
  inertia[1][1] = 1;
  inertia[2][2] = 1;

  fluid_forces[0] = 0;
  fluid_forces[1] = 0;

  velocity[0] = 0;
  velocity[1] = 0;

  position[0] = 0;
  position[1] = 0;

  fluid_torque[0] = 0;
  fluid_torque[1] = 0;
  fluid_torque[2] = 0;

  orientation[0] = 0;
  orientation[1] = 0;
  orientation[2] = 0;

  omega[0] = 0;
  omega[1] = 0;
  omega[2] = 0;

  if (dim == 3)
    {
      fluid_forces[2] = 0;
      velocity[2]     = 0;
      position[2]     = 0;
    }

  // Fill the vectors with default value
  previous_fluid_forces = fluid_forces;
  previous_fluid_torque = fluid_torque;
  velocity_iter         = velocity;

  omega_iter           = omega;
  omega_impulsion      = 0;
  omega_impulsion_iter = 0;
  impulsion            = 0;
  impulsion_iter       = 0;
  contact_impulsion    = 0;

  previous_positions.resize(3);
  previous_velocity.resize(3);
  previous_orientation.resize(3);
  previous_omega.resize(3);

  for (unsigned int i = 0; i < 3; ++i)
    {
      previous_positions[i]   = position;
      previous_velocity[i]    = velocity;
      previous_orientation[i] = orientation;
      previous_omega[i]       = omega;
    }
  residual_velocity = DBL_MAX;
  residual_omega    = DBL_MAX;
}

template <int dim>
void
IBParticle<dim>::initialise_last()
{
  // initialise all the variables associated to an immersed boundary particle
  previous_fluid_forces = fluid_forces;
  velocity_iter         = velocity;
  impulsion_iter        = impulsion;
  omega_iter            = omega;
  omega_impulsion_iter  = omega_impulsion;

  for (unsigned int i = 0; i < 3; ++i)
    {
      previous_positions[i]   = position;
      previous_velocity[i]    = velocity;
      previous_orientation[i] = orientation;
      previous_omega[i]       = omega;
    }
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
  properties[PropertiesIndex::m]  = std::make_pair("Mass", 1);
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
  properties[1]  = radius * 2.0;
  properties[2]  = velocity[0];
  properties[3]  = velocity[1];
  properties[5]  = fluid_forces[0];
  properties[6]  = fluid_forces[1];
  properties[8]  = omega[0];
  properties[9]  = omega[1];
  properties[10] = omega[2];
  properties[11] = mass;
  properties[12] = fluid_torque[0];
  properties[13] = fluid_torque[1];
  properties[14] = fluid_torque[2];
  if (dim == 2)
    {
      properties[4] = 0;
      properties[7] = 0;
    }
  if (dim == 3)
    {
      properties[4] = velocity[2];
      properties[7] = fluid_forces[2];
    }

  return properties;
}

template <int dim>
unsigned int
IBParticle<dim>::get_number_properties()
{
  return PropertiesIndex::n_properties;
}

template <int dim>
double
IBParticle<dim>::get_levelset(const Point<dim> &p)
{
  return shape->value(p);
}

template <int dim>
void
IBParticle<dim>::closest_surface_point(const Point<dim> &p,
                                       Point<dim> &      closest_point)
{
  Tensor<1, dim> actual_gradient       = shape->gradient(p);
  double         distance_from_surface = shape->value(p);
  closest_point =
    p - (actual_gradient / actual_gradient.norm()) * distance_from_surface;
}

template <int dim>
bool
IBParticle<dim>::is_inside_crown(const Point<dim> &evaluation_pt,
                                 const double      outer_radius,
                                 const double      inside_radius)
{
  bool is_inside;

  const double radius = shape->effective_radius;

  double distance              = shape->value(evaluation_pt);
  bool   is_inside_outer_ring  = distance <= radius * (outer_radius - 1);
  bool   is_outside_inner_ring = distance >= radius * (inside_radius - 1);

  is_inside = is_inside_outer_ring && is_outside_inner_ring;
  return is_inside;
}

template <int dim>
void
IBParticle<dim>::set_position(const Point<dim> position)
{
  this->position = position;
  this->shape->set_position(this->position);
}

template <int dim>
void
IBParticle<dim>::set_position(const double       position_component,
                              const unsigned int component)
{
  this->position[component] = position_component;
  this->shape->set_position(this->position);
}

template <int dim>
void
IBParticle<dim>::set_orientation(const Tensor<1, 3> orientation)
{
  this->orientation = orientation;
  this->shape->set_orientation(this->orientation);
}

template class IBParticle<2>;
template class IBParticle<3>;
