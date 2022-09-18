#include <core/ib_particle.h>
#include <core/shape.h>

#include <cfloat>


template <int dim>
void
IBParticle<dim>::initialize_all()
{
  // initialize all the variables associated to an immersed boundary particle.
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

  f_position    = std::make_shared<Functions::ParsedFunction<dim>>(dim);
  f_velocity    = std::make_shared<Functions::ParsedFunction<dim>>(dim);
  f_omega       = std::make_shared<Functions::ParsedFunction<dim>>(3);
  f_orientation = std::make_shared<Functions::ParsedFunction<dim>>(3);
}

template <int dim>
void
IBParticle<dim>::initialize_previous_solution()
{
  // initialize all the variables associated to an immersed boundary particle
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
void
IBParticle<dim>::initialize_shape(const std::string         type,
                                  const std::vector<double> shape_arguments)
{
  if (type == "sphere")
    shape =
      std::make_shared<Sphere<dim>>(shape_arguments[0], position, orientation);
  else if (type == "rectangle")
    {
      Tensor<1, dim> half_lengths;
      for (unsigned int i = 0; i < dim; ++i)
        {
          half_lengths[i] = shape_arguments[i];
        }
      shape =
        std::make_shared<Rectangle<dim>>(half_lengths, position, orientation);
    }
  else if (type == "ellipsoid")
    {
      Tensor<1, dim> radii;
      for (unsigned int i = 0; i < dim; ++i)
        {
          radii[i] = shape_arguments[i];
        }
      shape = std::make_shared<Ellipsoid<dim>>(radii, position, orientation);
    }
  else if (type == "torus")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<Torus<dim>>(shape_arguments[0],
                                             shape_arguments[1],
                                             position,
                                             orientation);
    }
  else if (type == "cone")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<Cone<dim>>(shape_arguments[0],
                                            shape_arguments[1],
                                            position,
                                            orientation);
    }
  else if (type == "cut hollow sphere")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<CutHollowSphere<dim>>(shape_arguments[0],
                                                       shape_arguments[1],
                                                       shape_arguments[2],
                                                       position,
                                                       orientation);
    }
  else if (type == "death star")
    {
      if constexpr (dim == 3)
        shape = std::make_shared<DeathStar<dim>>(shape_arguments[0],
                                                 shape_arguments[1],
                                                 shape_arguments[2],
                                                 position,
                                                 orientation);
    }
  else if (type == "rbf")
    {
      constexpr unsigned int numbers_per_nodes = dim + 3;
      unsigned int number_of_nodes = shape_arguments.size() / numbers_per_nodes;
      std::vector<double>         support_radius;
      std::vector<unsigned int>   basis_function;
      std::vector<double>         weight;
      std::vector<Tensor<1, dim>> nodes;
      support_radius.resize(number_of_nodes);
      basis_function.resize(number_of_nodes);
      weight.resize(number_of_nodes);
      nodes.resize(number_of_nodes);
      // The data is stored in this order: all weights, all radii, all basis
      // functions, all x positions, all y positions(, all z positions)
      for (unsigned int n_i = 0; n_i < number_of_nodes; n_i++)
        {
          weight[n_i]         = shape_arguments[0 * number_of_nodes + n_i];
          support_radius[n_i] = shape_arguments[1 * number_of_nodes + n_i];
          basis_function[n_i] =
            (unsigned int)round(shape_arguments[2 * number_of_nodes + n_i]);
          nodes[n_i][0] = shape_arguments[3 * number_of_nodes + n_i];
          nodes[n_i][1] = shape_arguments[4 * number_of_nodes + n_i];
          if constexpr (dim == 3)
            nodes[n_i][2] = shape_arguments[5 * number_of_nodes + n_i];
        }
      shape = std::make_shared<RBFShape<dim>>(
        support_radius, basis_function, weight, nodes, position, orientation);
    }
  else
    StandardExceptions::ExcNotImplemented();
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
IBParticle<dim>::is_inside_crown(const Point<dim> &evaluation_point,
                                 const double      outer_radius,
                                 const double      inside_radius)
{
  const double radius = shape->effective_radius;

  double distance              = shape->value(evaluation_point);
  bool   is_inside_outer_ring  = distance <= radius * (outer_radius - 1);
  bool   is_outside_inner_ring = distance >= radius * (inside_radius - 1);

  return is_inside_outer_ring && is_outside_inner_ring;
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
