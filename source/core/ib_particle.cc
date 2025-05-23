// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/ib_particle.h>
#include <core/shape_parsing.h>

// Deal.II includes
#include <deal.II/lac/lapack_full_matrix.h>

#include <cfloat>

template <int dim>
void
IBParticle<dim>::initialize_all()
{
  // initialize all the variables associated to an immersed boundary particle.
  radius      = 1.;
  particle_id = 0.;

  inertia[0][0] = 1.;
  inertia[0][1] = 0.;
  inertia[0][2] = 0.;
  inertia[1][0] = 0.;
  inertia[1][1] = 1.;
  inertia[1][2] = 0.;
  inertia[2][0] = 0.;
  inertia[2][1] = 0.;
  inertia[2][2] = 1.;


  fluid_forces[0] = 0.;
  fluid_forces[1] = 0.;

  velocity[0] = 0.;
  velocity[1] = 0.;

  position[0] = 0.;
  position[1] = 0.;

  fluid_torque[0] = 0.;
  fluid_torque[1] = 0.;
  fluid_torque[2] = 0.;

  orientation[0] = 0.;
  orientation[1] = 0.;
  orientation[2] = 0.;

  omega[0] = 0.;
  omega[1] = 0.;
  omega[2] = 0.;

  if (dim == 3)
    {
      fluid_forces[2] = 0.;
      velocity[2]     = 0.;
      position[2]     = 0.;
    }

  // Fill the vectors with default value
  previous_fluid_forces = fluid_forces;
  previous_fluid_torque = fluid_torque;
  velocity_iter         = velocity;

  omega_iter           = omega;
  omega_impulsion      = 0.;
  omega_impulsion_iter = 0.;
  impulsion            = 0.;
  impulsion_iter       = 0.;
  contact_impulsion    = 0.;

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

  f_position    = std::make_shared<Functions::ParsedFunction<dim>>(dim);
  f_velocity    = std::make_shared<Functions::ParsedFunction<dim>>(dim);
  f_omega       = std::make_shared<Functions::ParsedFunction<dim>>(3);
  f_orientation = std::make_shared<Functions::ParsedFunction<dim>>(3);

  rotation_matrix[0][0] = 1.;
  rotation_matrix[0][1] = 0.;
  rotation_matrix[0][2] = 0.;
  rotation_matrix[1][0] = 0.;
  rotation_matrix[1][1] = 1.;
  rotation_matrix[1][2] = 0.;
  rotation_matrix[2][0] = 0.;
  rotation_matrix[2][1] = 0.;
  rotation_matrix[2][2] = 1.;
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
  std::vector<std::pair<std::string, int>> particle_properties(
    DEM::CFDDEMProperties::PropertiesIndex::n_properties);
  particle_properties =
    DEM::ParticleProperties<dim, DEM::CFDDEMProperties::PropertiesIndex>::
      get_properties_name();

  // Rename particle type to id to match IB pattern
  particle_properties[DEM::CFDDEMProperties::type] = std::make_pair("id", 1);

  return particle_properties;
}

template <int dim>
std::vector<double>
IBParticle<dim>::get_properties()
{
  std::vector<double> properties(get_number_properties());
  properties[DEM::CFDDEMProperties::PropertiesIndex::type]    = particle_id;
  properties[DEM::CFDDEMProperties::PropertiesIndex::dp]      = radius * 2.0;
  properties[DEM::CFDDEMProperties::PropertiesIndex::v_x]     = velocity[0];
  properties[DEM::CFDDEMProperties::PropertiesIndex::v_y]     = velocity[1];
  properties[DEM::CFDDEMProperties::PropertiesIndex::omega_x] = omega[0];
  properties[DEM::CFDDEMProperties::PropertiesIndex::omega_y] = omega[1];
  properties[DEM::CFDDEMProperties::PropertiesIndex::omega_z] = omega[2];
  properties[DEM::CFDDEMProperties::PropertiesIndex::fem_force_x] =
    fluid_forces[0];
  properties[DEM::CFDDEMProperties::PropertiesIndex::fem_force_y] =
    fluid_forces[1];
  properties[DEM::CFDDEMProperties::PropertiesIndex::fem_torque_x] =
    fluid_torque[0];
  properties[DEM::CFDDEMProperties::PropertiesIndex::fem_torque_y] =
    fluid_torque[1];
  properties[DEM::CFDDEMProperties::PropertiesIndex::fem_torque_z] =
    fluid_torque[2];
  properties[DEM::CFDDEMProperties::PropertiesIndex::mass] = mass;

  if (dim == 2)
    {
      properties[DEM::CFDDEMProperties::PropertiesIndex::v_z]         = 0.;
      properties[DEM::CFDDEMProperties::PropertiesIndex::fem_force_z] = 0.;
    }
  if (dim == 3)
    {
      properties[DEM::CFDDEMProperties::PropertiesIndex::v_z] = velocity[2];
      properties[DEM::CFDDEMProperties::PropertiesIndex::fem_force_z] =
        fluid_forces[2];
    }

  return properties;
}

template <int dim>
void
IBParticle<dim>::clear_shape_cache()
{
  this->shape->clear_cache();
}

template <int dim>
void
IBParticle<dim>::initialize_shape(const std::string         type,
                                  const std::vector<double> shape_arguments)
{
  shape = ShapeGenerator::initialize_shape_from_vector(type,
                                                       shape_arguments,
                                                       position,
                                                       orientation);
}

template <int dim>
void
IBParticle<dim>::initialize_shape(const std::string type,
                                  const std::string raw_arguments)
{
  shape = ShapeGenerator::initialize_shape(type,
                                           raw_arguments,
                                           position,
                                           orientation);
}

template <int dim>
void
IBParticle<dim>::closest_surface_point(
  const Point<dim>                                     &p,
  Point<dim>                                           &closest_point,
  const typename DoFHandler<dim>::active_cell_iterator &cell_guess)
{
  shape->closest_surface_point(p, closest_point, cell_guess);
}

template <int dim>
void
IBParticle<dim>::closest_surface_point(const Point<dim> &p,
                                       Point<dim>       &closest_point)
{
  shape->closest_surface_point(p, closest_point);
}

template <int dim>
bool
IBParticle<dim>::is_inside_crown(
  const Point<dim>                                     &evaluation_point,
  const double                                          outer_radius,
  const double                                          inside_radius,
  const bool                                            absolute_distance,
  const typename DoFHandler<dim>::active_cell_iterator &cell_guess)
{
  const double radius = shape->effective_radius;

  double distance = shape->value_with_cell_guess(evaluation_point, cell_guess);
  bool   is_inside_outer_ring;
  bool   is_outside_inner_ring;
  if (absolute_distance)
    {
      is_inside_outer_ring  = distance <= outer_radius;
      is_outside_inner_ring = distance >= inside_radius;
    }
  else
    {
      is_inside_outer_ring  = distance <= radius * (outer_radius - 1);
      is_outside_inner_ring = distance >= radius * (inside_radius - 1);
    }

  return is_inside_outer_ring && is_outside_inner_ring;
}

template <int dim>
bool
IBParticle<dim>::is_inside_crown(const Point<dim> &evaluation_point,
                                 const double      outer_radius,
                                 const double      inside_radius,
                                 const bool        absolute_distance)
{
  const double radius = shape->effective_radius;

  double distance = shape->value(evaluation_point);
  bool   is_inside_outer_ring;
  bool   is_outside_inner_ring;
  if (absolute_distance)
    {
      is_inside_outer_ring  = distance <= outer_radius;
      is_outside_inner_ring = distance >= inside_radius;
    }
  else
    {
      is_inside_outer_ring  = distance <= radius * (outer_radius - 1);
      is_outside_inner_ring = distance >= radius * (inside_radius - 1);
    }

  return is_inside_outer_ring && is_outside_inner_ring;
}

template <int dim>
void
IBParticle<dim>::set_orientation(const Tensor<1, 3> new_orientation)
{
  this->orientation = new_orientation;
  this->shape->set_orientation(new_orientation);
  this->rotation_matrix = this->shape->get_rotation_matrix();
}

template <int dim>
void
IBParticle<dim>::update_precalculations(DoFHandler<dim> &updated_dof_handler,
                                        const bool mesh_based_precalculations)
{
  if (integrate_motion || velocity.norm() > 1e-16 || omega.norm() > 1e-16)
    {
      this->load_data_from_file();
    }
  if (typeid(*shape) == typeid(RBFShape<dim>))
    {
      std::static_pointer_cast<RBFShape<dim>>(shape)->update_precalculations(
        updated_dof_handler, mesh_based_precalculations);
    }
  else if (typeid(*shape) == typeid(CompositeShape<dim>))
    {
      std::static_pointer_cast<CompositeShape<dim>>(shape)
        ->update_precalculations(updated_dof_handler,
                                 mesh_based_precalculations);
    }
}

template <int dim>
void
IBParticle<dim>::remove_superfluous_data(DoFHandler<dim> &updated_dof_handler,
                                         const bool mesh_based_precalculations)
{
  if (typeid(*shape) == typeid(RBFShape<dim>))
    std::static_pointer_cast<RBFShape<dim>>(shape)->remove_superfluous_data(
      updated_dof_handler, mesh_based_precalculations);
  else if (typeid(*shape) == typeid(CompositeShape<dim>))
    std::static_pointer_cast<CompositeShape<dim>>(shape)
      ->remove_superfluous_data(updated_dof_handler,
                                mesh_based_precalculations);
}

template <int dim>
void
IBParticle<dim>::load_data_from_file()
{
  if (typeid(*shape) == typeid(RBFShape<dim>))
    std::static_pointer_cast<RBFShape<dim>>(shape)->load_data_from_file();
  else if (typeid(*shape) == typeid(CompositeShape<dim>))
    std::static_pointer_cast<CompositeShape<dim>>(shape)->load_data_from_file();
}

template class IBParticle<2>;
template class IBParticle<3>;
