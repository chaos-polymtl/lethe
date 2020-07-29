/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <dem/pp_contact_force.h>

// Updates the contact information (contact_info) based on the new information
// of particles pair in the current time step
template <int dim>
void
PPContactForce<dim>::update_contact_information(
  pp_contact_info_struct<dim> &  contact_info,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const Point<dim> &             particle_one_location,
  const Point<dim> &             particle_two_location,
  const double &                 dt)
{
  // Calculation of the contact vector (vector from particle one to particle two
  auto contact_vector = particle_two_location - particle_one_location;

  // Using contact_vector, the contact normal vector is obtained
  auto normal_unit_vector = contact_vector / contact_vector.norm();

  // Defining velocities and angular velocities of particles one and
  // two as vectors
  Tensor<1, dim> particle_one_velocity, particle_two_velocity,
    particle_one_omega, particle_two_omega;

  // Defining relative contact velocity
  Tensor<1, dim> contact_relative_velocity;

  // Finding velocities and angular velocities of particles
  particle_one_velocity[0] = particle_one_properties[DEM::PropertiesIndex::v_x];
  particle_one_velocity[1] = particle_one_properties[DEM::PropertiesIndex::v_y];

  particle_two_velocity[0] = particle_two_properties[DEM::PropertiesIndex::v_x];
  particle_two_velocity[1] = particle_two_properties[DEM::PropertiesIndex::v_y];

  particle_one_omega[0] =
    particle_one_properties[DEM::PropertiesIndex::omega_x];
  particle_one_omega[1] =
    particle_one_properties[DEM::PropertiesIndex::omega_y];

  particle_two_omega[0] =
    particle_two_properties[DEM::PropertiesIndex::omega_x];
  particle_two_omega[1] =
    particle_two_properties[DEM::PropertiesIndex::omega_y];

  if (dim == 3)
    {
      particle_one_velocity[2] =
        particle_one_properties[DEM::PropertiesIndex::v_z];
      particle_two_velocity[2] =
        particle_two_properties[DEM::PropertiesIndex::v_z];
      particle_one_omega[2] =
        particle_one_properties[DEM::PropertiesIndex::omega_z];
      particle_two_omega[2] =
        particle_two_properties[DEM::PropertiesIndex::omega_z];

      // Calculation of contact relative velocity
      contact_relative_velocity =
        (particle_one_velocity - particle_two_velocity) +
        (cross_product_3d(0.5 *
                            (particle_one_properties[DEM::PropertiesIndex::dp] *
                               particle_one_omega +
                             particle_two_properties[DEM::PropertiesIndex::dp] *
                               particle_two_omega),
                          normal_unit_vector));
    }
  else
    {
      // if dim == 2
      contact_relative_velocity = particle_one_velocity - particle_two_velocity;
    }

  // Calculation of normal relative velocity. Note that in the
  // following line the product acts as inner product since both
  // sides are vectors, while in the second line the product is
  // scalar and vector product
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_unit_vector;
  Tensor<1, dim> normal_relative_velocity =
    normal_relative_velocity_value * normal_unit_vector;

  // Calculation of tangential relative velocity
  Tensor<1, dim> tangential_relative_velocity =
    contact_relative_velocity - normal_relative_velocity;

  // Calculation of new tangential_overlap, since this value is
  // history-dependent it needs the value at previous time-step
  // This variable is the main reason that we have iteration over
  // two different vectors (pairs_in_contact and
  // contact_pair_candidates): tangential_overlap of the particles
  // which were already in contact (pairs_in_contact) needs to
  // modified using its history, while the tangential_overlaps of
  // new particles are equal to zero
  Tensor<1, dim> last_step_tangential_overlap = contact_info.tangential_overlap;
  Tensor<1, dim> tangential_overlap =
    last_step_tangential_overlap -
    (last_step_tangential_overlap * normal_unit_vector) * normal_unit_vector;

  // Adding a small value to the tangential_overlap_norm to avoid
  // 0/0 occurance
  double tangential_overlap_norm = tangential_overlap.norm() + DBL_MIN;

  Tensor<1, dim> modified_tangential_overlap =
    (last_step_tangential_overlap.norm() / tangential_overlap_norm) *
      tangential_overlap +
    contact_info.tangential_relative_velocity * dt;

  // std::cout << contact_info.tangential_relative_velocity << std::endl;
  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.normal_unit_vector           = normal_unit_vector;
  contact_info.tangential_overlap           = modified_tangential_overlap;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
}

// This function is used to apply calculated forces and torques on the particle
// pair
template <int dim>
void
PPContactForce<dim>::apply_force_and_torque_real(
  ArrayView<double> &particle_one_properties,
  ArrayView<double> &particle_two_properties,
  const std::
    tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
      &forces_and_torques)
{
  // Getting the values from the forces_and_torques tuple, which are: 1, normal
  // force, 2, tangential force, 3, tangential torque and 4, rolling resistance
  // torque
  Tensor<1, dim> normal_force              = std::get<0>(forces_and_torques);
  Tensor<1, dim> tangential_force          = std::get<1>(forces_and_torques);
  Tensor<1, dim> tangential_torque         = std::get<2>(forces_and_torques);
  Tensor<1, dim> rolling_resistance_torque = std::get<3>(forces_and_torques);

  // Calculation of total force
  Tensor<1, dim> total_force = normal_force + tangential_force;

  // Updating the force of particles in the particle handler
  for (int d = 0; d < dim; ++d)
    {
      particle_one_properties[DEM::PropertiesIndex::force_x + d] =
        particle_one_properties[DEM::PropertiesIndex::force_x + d] -
        total_force[d];
      particle_two_properties[DEM::PropertiesIndex::force_x + d] =
        particle_two_properties[DEM::PropertiesIndex::force_x + d] +
        total_force[d];
    }

  // Updating the torque acting on particles
  for (int d = 0; d < dim; ++d)
    {
      particle_one_properties[DEM::PropertiesIndex::M_x + d] =
        particle_one_properties[DEM::PropertiesIndex::M_x + d] -
        tangential_torque[d] + rolling_resistance_torque[d];
      particle_two_properties[DEM::PropertiesIndex::M_x + d] =
        particle_two_properties[DEM::PropertiesIndex::M_x + d] -
        tangential_torque[d] - rolling_resistance_torque[d];
    }
}

// This function is used to apply calculated forces and torques on the particle
// pair
template <int dim>
void
PPContactForce<dim>::apply_force_and_torque_ghost(
  ArrayView<double> &particle_one_properties,
  const std::
    tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
      &forces_and_torques)
{
  // Getting the values from the forces_and_torques tuple, which are: 1, normal
  // force, 2, tangential force, 3, tangential torque and 4, rolling resistance
  // torque
  Tensor<1, dim> normal_force              = std::get<0>(forces_and_torques);
  Tensor<1, dim> tangential_force          = std::get<1>(forces_and_torques);
  Tensor<1, dim> tangential_torque         = std::get<2>(forces_and_torques);
  Tensor<1, dim> rolling_resistance_torque = std::get<3>(forces_and_torques);

  // Calculation of total force
  Tensor<1, dim> total_force = normal_force + tangential_force;

  // Updating the force of particles in the particle handler
  for (int d = 0; d < dim; ++d)
    {
      particle_one_properties[DEM::PropertiesIndex::force_x + d] =
        particle_one_properties[DEM::PropertiesIndex::force_x + d] -
        total_force[d];
    }

  // Updating the torque acting on particles
  for (int d = 0; d < dim; ++d)
    {
      particle_one_properties[DEM::PropertiesIndex::M_x + d] =
        particle_one_properties[DEM::PropertiesIndex::M_x + d] -
        tangential_torque[d] + rolling_resistance_torque[d];
    }
}

template class PPContactForce<2>;
template class PPContactForce<3>;
