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

#include <dem/pw_contact_force.h>

// Updates the contact information (contact_info) based on the new information
// of particles pair in the current time step
template <int dim>
void
PWContactForce<dim>::update_contact_information(
  pw_contact_info_struct<dim> &  contact_info,
  const ArrayView<const double> &particle_properties,
  const double &                 dt)
{
  auto normal_vector = contact_info.normal_vector;

  // Using velocity and angular velocity of particle as
  // local vectors
  Tensor<1, dim> particle_velocity;
  particle_velocity[0] = particle_properties[DEM::PropertiesIndex::v_x];
  particle_velocity[1] = particle_properties[DEM::PropertiesIndex::v_y];
  if (dim == 3)
    {
      particle_velocity[2] = particle_properties[DEM::PropertiesIndex::v_z];
    }

  Tensor<1, dim> particle_omega;
  particle_omega[0] = particle_properties[DEM::PropertiesIndex::omega_x];
  particle_omega[1] = particle_properties[DEM::PropertiesIndex::omega_y];
  if (dim == 3)
    {
      particle_omega[2] = particle_properties[DEM::PropertiesIndex::omega_z];
    }

  // Defining relative contact velocity
  Tensor<1, dim> contact_relative_velocity;
  if (dim == 3)
    {
      contact_relative_velocity =
        particle_velocity +
        cross_product_3d((((particle_properties[DEM::PropertiesIndex::dp]) /
                           2) *
                          particle_omega),
                         normal_vector);
    }

  if (dim == 2)
    {
      contact_relative_velocity = particle_velocity;
    }

  // Calculation of normal relative velocity
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_vector;
  Tensor<1, dim> normal_relative_velocity =
    normal_relative_velocity_value * normal_vector;

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
    (last_step_tangential_overlap * normal_vector) * normal_vector;

  // Adding a small value to the tangential_overlap_norm to avoid
  // 0/0 occurance
  double tangential_overlap_norm = tangential_overlap.norm() + DBL_MIN;

  Tensor<1, dim> modified_tangential_overlap =
    (last_step_tangential_overlap.norm() / tangential_overlap_norm) *
      tangential_overlap +
    contact_info.tangential_relative_velocity * dt;

  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.tangential_overlap           = modified_tangential_overlap;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
}

// This function is used to apply calculated forces and torques on the particle
// pair
template <int dim>
void
PWContactForce<dim>::apply_force_and_torque(
  ArrayView<double> &particle_properties,
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
      particle_properties[DEM::PropertiesIndex::force_x + d] =
        particle_properties[DEM::PropertiesIndex::force_x + d] + total_force[d];
    }

  // Updating the torque acting on particles
  for (int d = 0; d < dim; ++d)
    {
      particle_properties[DEM::PropertiesIndex::M_x + d] =
        particle_properties[DEM::PropertiesIndex::M_x + d] +
        tangential_torque[d] + rolling_resistance_torque[d];
    }
}

template <int dim>
Tensor<1, dim> PWContactForce<dim>::find_projection(Tensor<1, dim> vector_a,
                                                    Tensor<1, dim> vector_b)
{
  Tensor<1, dim> vector_c;
  vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

  return vector_c;
}

template class PWContactForce<2>;
template class PWContactForce<3>;
