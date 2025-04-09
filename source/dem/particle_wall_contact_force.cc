// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/particle_wall_contact_force.h>

// Updates the contact information (contact_info) based on the new information
// of particles pair in the current time step
template <int dim, typename PropertiesIndex>
void
ParticleWallContactForce<dim, PropertiesIndex>::update_contact_information(
  particle_wall_contact_info<dim> &contact_info,
  const Point<3>                  &particle_position,
  const ArrayView<const double>   &particle_properties,
  const double                     dt)
{
  // i is the particle, j is the wall.
  // we need to put a minus sign infront of the normal_vector to respect the
  // convention (i -> j)
  auto               normal_vector = -contact_info.normal_vector;
  const unsigned int boundary_id   = contact_info.boundary_id;

  // Using velocity and angular velocity of particle as
  // local vectors
  Tensor<1, 3> particle_velocity;
  particle_velocity[0] = particle_properties[PropertiesIndex::v_x];
  particle_velocity[1] = particle_properties[PropertiesIndex::v_y];
  particle_velocity[2] = particle_properties[PropertiesIndex::v_z];

  Tensor<1, 3> particle_angular_velocity;
  particle_angular_velocity[0] = particle_properties[PropertiesIndex::omega_x];
  particle_angular_velocity[1] = particle_properties[PropertiesIndex::omega_y];
  particle_angular_velocity[2] = particle_properties[PropertiesIndex::omega_z];

  // Calculate approximation of the contact point using the normal vector
  Point<3> contact_point =
    particle_position +
    0.5 * particle_properties[PropertiesIndex::dp] * normal_vector;

  // Get vector pointing from the contact point to the origin of the rotation
  // axis
  Tensor<1, 3> vector_to_rotating_axis =
    contact_point - this->point_on_rotation_vector[boundary_id];

  // Remove the rotating axis component of that vector
  vector_to_rotating_axis =
    vector_to_rotating_axis -
    (vector_to_rotating_axis * this->boundary_rotational_vector[boundary_id]) *
      this->boundary_rotational_vector[boundary_id];

  // Tensor<1,3> vector_to_rotation_axis = this->boundary_rotational_speed_map

  // Defining relative contact velocity using the convention
  // v_ij = v_j - v_i
  Tensor<1, 3> contact_relative_velocity =
    this->boundary_translational_velocity_map[boundary_id] - particle_velocity +
    cross_product_3d((-0.5 * particle_properties[PropertiesIndex::dp] *
                      particle_angular_velocity),
                     normal_vector) +
    cross_product_3d(this->boundary_rotational_speed_map[boundary_id] *
                       this->boundary_rotational_vector[boundary_id],
                     vector_to_rotating_axis);

  // Calculation of normal relative velocity
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_vector;
  Tensor<1, 3> normal_relative_velocity =
    normal_relative_velocity_value * normal_vector;

  // Calculation of tangential relative velocity
  Tensor<1, 3> tangential_relative_velocity =
    contact_relative_velocity - normal_relative_velocity;

  // Calculation of new tangential_displacement, since this value is
  // history-dependent it needs the value at previous time-step
  // This variable is the main reason that we have iteration over
  // two different vectors (pairs_in_contact and
  // contact_pair_candidates): tangential_displacement of the particles
  // which were already in contact (pairs_in_contact) needs to be
  // modified using its history, while the tangential_displacements of
  // new particles are equal to zero
  Tensor<1, 3> modified_tangential_displacement =
    contact_info.tangential_displacement + tangential_relative_velocity * dt;

  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.tangential_displacement           = modified_tangential_displacement;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
}

template <int dim, typename PropertiesIndex>
void
ParticleWallContactForce<dim, PropertiesIndex>::
  update_particle_floating_wall_contact_information(
    particle_wall_contact_info<dim> &contact_info,
    const ArrayView<const double>   &particle_properties,
    const double                     dt,
    const Tensor<1, 3>              &cut_cell_translational_velocity,
    const Tensor<1, 3>              &cut_cell_rotational_velocity,
    const double                     center_of_rotation_particle_distance)
{
  // i is the particle, j is the wall.
  // we need to put a minus sign infront of the normal_vector to respect the
  // convention (i -> j)
  const Tensor<1, 3> normal_vector = -contact_info.normal_vector;

  // Using velocity and angular velocity of particle as
  // local vectors
  Tensor<1, 3> particle_velocity;
  particle_velocity[0] = particle_properties[PropertiesIndex::v_x];
  particle_velocity[1] = particle_properties[PropertiesIndex::v_y];
  particle_velocity[2] = particle_properties[PropertiesIndex::v_z];


  Tensor<1, 3> particle_angular_velocity;
  particle_angular_velocity[0] = particle_properties[PropertiesIndex::omega_x];
  particle_angular_velocity[1] = particle_properties[PropertiesIndex::omega_y];
  particle_angular_velocity[2] = particle_properties[PropertiesIndex::omega_z];

  // Defining relative contact velocity
  // v_ij = v_j - v_i
  Tensor<1, 3> contact_relative_velocity =
    cut_cell_translational_velocity - particle_velocity +
    cross_product_3d((center_of_rotation_particle_distance *
                        cut_cell_rotational_velocity -
                      0.5 * particle_properties[PropertiesIndex::dp] *
                        particle_angular_velocity),
                     normal_vector);

  // Calculation of normal relative velocity
  double normal_relative_velocity_value =
    contact_relative_velocity * normal_vector;
  Tensor<1, 3> normal_relative_velocity =
    normal_relative_velocity_value * normal_vector;

  // Calculation of tangential relative velocity
  Tensor<1, 3> tangential_relative_velocity =
    contact_relative_velocity - normal_relative_velocity;

  // Calculate the new tangential_displacement, since this value is
  // history-dependent, it needs the value at previous time-step
  Tensor<1, 3> modified_tangential_displacement =
    contact_info.tangential_displacement + tangential_relative_velocity * dt;

  // Updating the contact_info container based on the new calculated values
  contact_info.normal_relative_velocity     = normal_relative_velocity_value;
  contact_info.tangential_displacement           = modified_tangential_displacement;
  contact_info.tangential_relative_velocity = tangential_relative_velocity;
}

template <int dim, typename PropertiesIndex>
void
ParticleWallContactForce<dim, PropertiesIndex>::
  calculate_force_and_torque_on_boundary(const unsigned int boundary_id,
                                         Tensor<1, 3>       add_force,
                                         const Point<3>     point_contact)
{
  if (calculate_force_torque_on_boundary)
    {
      mpi_correction_over_calculation_of_forces_and_torques();

      force_on_walls[boundary_id] = force_on_walls[boundary_id] - add_force;

      torque_on_walls[boundary_id] =
        torque_on_walls[boundary_id] -
        cross_product_3d(point_contact - center_mass_container, add_force);
    }
}

template <int dim, typename PropertiesIndex>
std::map<unsigned int, Tensor<1, 3>>
ParticleWallContactForce<dim, PropertiesIndex>::initialize()
{
  std::map<unsigned int, Tensor<1, 3>> map;
  for (const auto &it : boundary_index)
    {
      map[it] = 0;
    }
  return map;
}

template <int dim, typename PropertiesIndex>
void
ParticleWallContactForce<dim, PropertiesIndex>::
  mpi_correction_over_calculation_of_forces_and_torques()
{
  for (const auto &it : boundary_index)
    {
      force_on_walls[it] =
        Utilities::MPI::sum(force_on_walls[it], MPI_COMM_WORLD);
      torque_on_walls[it] =
        Utilities::MPI::sum(torque_on_walls[it], MPI_COMM_WORLD);
    }
}


template class ParticleWallContactForce<2, DEM::DEMProperties::PropertiesIndex>;
template class ParticleWallContactForce<2,
                                        DEM::CFDDEMProperties::PropertiesIndex>;
template class ParticleWallContactForce<2,
                                        DEM::DEMMPProperties::PropertiesIndex>;
template class ParticleWallContactForce<3, DEM::DEMProperties::PropertiesIndex>;
template class ParticleWallContactForce<3,
                                        DEM::CFDDEMProperties::PropertiesIndex>;
template class ParticleWallContactForce<3,
                                        DEM::DEMMPProperties::PropertiesIndex>;
