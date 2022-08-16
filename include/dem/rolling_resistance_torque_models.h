/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 */

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/particles/particle_handler.h>

#ifndef rolling_resistance_torque_models_h
#  define rolling_resistance_torque_models_h

enum class RollingResistanceTorqueModel
{
  no_rolling_resistance,
  constant_rolling_resistance,
  viscous_rolling_resistance
};


/**
 * @brief
 */
inline Tensor<1, 3>
no_rolling_resistance_torque(
  const double /*effective_r*/,
  const ArrayView<const double> & /*particle_one_properties*/,
  const ArrayView<const double> & /*particle_two_properties*/,
  const double /*effective_rolling_friction_coefficient*/,
  const double /*normal_force_norm*/,
  const Tensor<1, 3> & /*normal_contact_vector*/)

{
  // No rolling resistance torque model. When this model is use, the rolling
  // friction is zero.
  Tensor<1, 3> rolling_resistance_torque({0.0, 0.0, 0.0});

  return rolling_resistance_torque;
}

/**
 * @brief
 */
// Constant resistance model
// M_r = - mu_r * R_eff * |F_n| * omega_hat
// omega_hat = (omega_i - omega_j) / (|oemaga_i - omega_j|)
inline Tensor<1, 3>
constant_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   normal_force_norm,
  const Tensor<1, 3> & /*normal_contact_vector*/)

{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity,
    omega_ij, omega_ij_direction;
  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  // Calculation of rolling resistance torque
  return (-effective_rolling_friction_coefficient * effective_r *
          normal_force_norm * omega_ij_direction);
}

/**
 * @brief
 */
// Viscous resistance model
// M_r = - mu_r * R_eff * |F_n| * |V_omega| * omega_hat
// omega_hat = (omega_i - omega_j) / (|oemaga_i - omega_j|)
// V_omega = omega_i × (R_i * n_ij) - omega_j × (R_j * n_ji)
inline Tensor<1, 3>
viscous_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   normal_force_norm,
  const Tensor<1, 3> &           normal_contact_vector)

{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity,
    omega_ij, omega_ij_direction;

  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex::omega_x + d];
    }

  omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
  omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  Tensor<1, 3> v_omega =
    cross_product_3d(particle_one_angular_velocity,
                     particle_one_properties[DEM::PropertiesIndex::dp] * 0.5 *
                       normal_contact_vector) -
    cross_product_3d(particle_two_angular_velocity,
                     particle_two_properties[DEM::PropertiesIndex::dp] * 0.5 *
                       -normal_contact_vector);

  // Calculation of rolling resistance torque
  return (-effective_rolling_friction_coefficient * effective_r *
          normal_force_norm * v_omega.norm() * omega_ij_direction);
}

#endif
