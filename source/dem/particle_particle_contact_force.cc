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

#include <dem/particle_particle_contact_force.h>

using namespace DEM;

// Updates the contact information (contact_info) based on the new
// information of particles pair in the current time step
template <int dim>
void
ParticleParticleContactForce<dim>::update_contact_information(
  particle_particle_contact_info_struct<dim> &contact_info,
  double &                                    normal_relative_velocity_value,
  Tensor<1, 3> &                              normal_unit_vector,
  const ArrayView<const double> &             particle_one_properties,
  const ArrayView<const double> &             particle_two_properties,
  const Point<3> &                            particle_one_location,
  const Point<3> &                            particle_two_location,
  const double                                dt)
{
  // Calculation of the contact vector (vector from particle one to particle two
  auto contact_vector = particle_two_location - particle_one_location;

  // Using contact_vector, the contact normal vector is obtained
  normal_unit_vector = contact_vector / contact_vector.norm();

  // Defining velocities and angular velocities of particles one and
  // two as vectors
  Tensor<1, 3> particle_one_velocity, particle_two_velocity, particle_one_omega,
    particle_two_omega;

  // Defining relative contact velocity
  Tensor<1, 3> contact_relative_velocity;

  // Finding velocities and angular velocities of particles
  particle_one_velocity[0] = particle_one_properties[PropertiesIndex::v_x];
  particle_one_velocity[1] = particle_one_properties[PropertiesIndex::v_y];

  particle_two_velocity[0] = particle_two_properties[PropertiesIndex::v_x];
  particle_two_velocity[1] = particle_two_properties[PropertiesIndex::v_y];

  particle_one_omega[0] = particle_one_properties[PropertiesIndex::omega_x];
  particle_one_omega[1] = particle_one_properties[PropertiesIndex::omega_y];

  particle_two_omega[0] = particle_two_properties[PropertiesIndex::omega_x];
  particle_two_omega[1] = particle_two_properties[PropertiesIndex::omega_y];

  particle_one_velocity[2] = particle_one_properties[PropertiesIndex::v_z];
  particle_two_velocity[2] = particle_two_properties[PropertiesIndex::v_z];
  particle_one_omega[2]    = particle_one_properties[PropertiesIndex::omega_z];
  particle_two_omega[2]    = particle_two_properties[PropertiesIndex::omega_z];

  // Calculation of contact relative velocity
  contact_relative_velocity =
    (particle_one_velocity - particle_two_velocity) +
    (cross_product_3d(
      0.5 * (particle_one_properties[PropertiesIndex::dp] * particle_one_omega +
             particle_two_properties[PropertiesIndex::dp] * particle_two_omega),
      normal_unit_vector));


  // Calculation of normal relative velocity. Note that in the
  // following line the product acts as inner product since both
  // sides are vectors, while in the second line the product is
  // scalar and vector product
  normal_relative_velocity_value =
    contact_relative_velocity * normal_unit_vector;

  // Calculation of tangential relative velocity
  contact_info.tangential_relative_velocity =
    contact_relative_velocity -
    (normal_relative_velocity_value * normal_unit_vector);

  // Calculation of new tangential_overlap, since this value is
  // history-dependent it needs the value at previous time-step
  // This variable is the main reason that we have iteration over
  // two different vectors (pairs_in_contact and
  // contact_pair_candidates): tangential_overlap of the particles
  // which were already in contact (pairs_in_contact) needs to
  // modified using its history, while the tangential_overlaps of
  // new particles are equal to zero
  contact_info.tangential_overlap +=
    contact_info.tangential_relative_velocity * dt;
}

template <int dim>
inline void
ParticleParticleContactForce<dim>::find_effective_radius_and_mass(
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties)
{
  effective_mass = (particle_one_properties[DEM::PropertiesIndex::mass] *
                    particle_two_properties[DEM::PropertiesIndex::mass]) /
                   (particle_one_properties[DEM::PropertiesIndex::mass] +
                    particle_two_properties[DEM::PropertiesIndex::mass]);
  effective_radius = (particle_one_properties[DEM::PropertiesIndex::dp] *
                      particle_two_properties[DEM::PropertiesIndex::dp]) /
                     (2 * (particle_one_properties[DEM::PropertiesIndex::dp] +
                           particle_two_properties[DEM::PropertiesIndex::dp]));
}

template class ParticleParticleContactForce<2>;
template class ParticleParticleContactForce<3>;
