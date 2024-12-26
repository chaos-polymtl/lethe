// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_rolling_resistance_torque_models_h
#define lethe_rolling_resistance_torque_models_h

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>

#include <deal.II/particles/particle_handler.h>

/**
 * @brief No rolling resistance torque model. No calculation is being done with
 * this model.
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
 * @brief Calculation of the constant rolling resistance torque model using the
 * information obtained from the fine search and physical properties.
 *
 * The model woks as follows:
 * M_r = - mu_r * R_eff * |F_n| * omega_hat
 * omega_hat = (omega_i - omega_j) / (|oemaga_i - omega_j|)
 * @param[in] effective_r Effective radius.
 * @param[in] particle_one_properties Properties of particle one in contact.
 * @param[in] particle_two_properties Properties of particle two in contact.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the nomal force.
 *
 */
template <DEM::SolverType solver_type>
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
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
    }

  Tensor<1, 3> omega_ij =
    particle_one_angular_velocity - particle_two_angular_velocity;
  Tensor<1, 3> omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  // Calculation of rolling resistance torque
  return (-effective_rolling_friction_coefficient * effective_r *
          normal_force_norm * omega_ij_direction);
}

/**
 * @brief Calculation of the viscous rolling resistance torque model using the
 * information obtained from the fine search and physical properties.
 *
 * The model woks as follows:
 * M_r = - mu_r * R_eff * |F_n| * |V_omega| * omega_hat
 * omega_hat = (omega_i - omega_j) / (|oemaga_i - omega_j|)
 * V_omega = omega_i × (R_i * n_ij) - omega_j × (R_j * n_ji)
 * @param[in] effective_r Effective radius.
 * @param[in] particle_one_properties Properties of particle one in contact.
 * @param[in] particle_two_properties Properties of particle two in contact.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the nomal force.
 * @param[in] normal_contact_vector Normal unit vector.
 */

template <DEM::SolverType solver_type>
inline Tensor<1, 3>
viscous_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   normal_force_norm,
  const Tensor<1, 3>            &normal_contact_vector)

{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
    }

  Tensor<1, 3> omega_ij =
    particle_one_angular_velocity - particle_two_angular_velocity;
  Tensor<1, 3> omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  Tensor<1, 3> v_omega =
    cross_product_3d(
      particle_one_angular_velocity,
      particle_one_properties[DEM::PropertiesIndex<solver_type>::dp] * 0.5 *
        normal_contact_vector) -
    cross_product_3d(
      particle_two_angular_velocity,
      particle_two_properties[DEM::PropertiesIndex<solver_type>::dp] * 0.5 *
        -normal_contact_vector);

  // Calculation of rolling resistance torque
  return (-effective_rolling_friction_coefficient * effective_r *
          normal_force_norm * v_omega.norm() * omega_ij_direction);
}

#endif
