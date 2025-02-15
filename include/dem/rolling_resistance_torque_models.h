// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_rolling_resistance_torque_models_h
#define lethe_rolling_resistance_torque_models_h

/**
 * @brief No rolling resistance torque model. No calculation is being done with
 * this model.
 */

inline Tensor<1, 3>
no_rolling_resistance_torque()
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
 * omega_hat = (omega_i - omega_j) / (|omega_i - omega_j|)
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] effective_r Effective radius.
 * @param[in] particle_one_properties Properties of particle one in contact.
 * @param[in] particle_two_properties Properties of particle two in contact.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the normal force.
 *
 */
template <typename PropertiesIndex>
inline Tensor<1, 3>
constant_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   normal_force_norm)

{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[PropertiesIndex::omega_x + d];
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
 * omega_hat = (omega_i - omega_j) / (|omega_i - omega_j|)
 * V_omega = omega_i × (R_i * n_ij) - omega_j × (R_j * n_ji)
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] effective_r Effective radius.
 * @param[in] particle_one_properties Properties of particle one in contact.
 * @param[in] particle_two_properties Properties of particle two in contact.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the normal force.
 * @param[in] normal_unit_vector Normal unit vector.
 */

template <typename PropertiesIndex>
inline Tensor<1, 3>
viscous_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   normal_force_norm,
  const Tensor<1, 3>            &normal_unit_vector)

{
  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  Tensor<1, 3> particle_one_angular_velocity, particle_two_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_one_angular_velocity[d] =
        particle_one_properties[PropertiesIndex::omega_x + d];
      particle_two_angular_velocity[d] =
        particle_two_properties[PropertiesIndex::omega_x + d];
    }

  Tensor<1, 3> omega_ij =
    particle_one_angular_velocity - particle_two_angular_velocity;
  Tensor<1, 3> omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

  Tensor<1, 3> v_omega =
    cross_product_3d(particle_one_angular_velocity,
                     particle_one_properties[PropertiesIndex::dp] * 0.5 *
                       normal_unit_vector) -
    cross_product_3d(particle_two_angular_velocity,
                     particle_two_properties[PropertiesIndex::dp] * 0.5 *
                       -normal_unit_vector);

  // Calculation of rolling resistance torque
  return (-effective_rolling_friction_coefficient * effective_r *
          normal_force_norm * v_omega.norm() * omega_ij_direction);
}

/**
 * @brief
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the problem is solved.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] effective_r Effective radius.
 * @param[in] particle_one_properties Properties of particle one in contact.
 * @param[in] particle_two_properties Properties of particle two in contact.
 * @param[in] effective_rolling_viscous_damping_coefficient
 * Effective rolling viscous damping
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] f_coefficient Model parameter for the viscous damping.
 * @param[in] normal_force_norm Norm of the normal force.
 * @param[in] dt DEM time step.
 * @param[in] normal_spring_constant normal contact stiffness constant.
 * @param[in] normal_unit_vector Normal unit vector between particles in
 * contact.
 * @param[in,out] cumulative_rolling_resistance_spring_torque Cumulative
 * rolling resistance torque for the EPSD rolling resistance model.
 */

template <int dim, typename PropertiesIndex>
inline Tensor<1, 3>
epsd_rolling_resistance_torque(
  const double                   effective_r,
  const ArrayView<const double> &particle_one_properties,
  const ArrayView<const double> &particle_two_properties,
  const double                   effective_rolling_friction_coefficient,
  const double                   effective_rolling_viscous_damping_coefficient,
  const double                   f_coefficient,
  const double                   normal_force_norm,
  const double                   dt,
  const double                   normal_spring_constant,
  const Tensor<1, 3>            &normal_unit_vector,
  Tensor<1, 3>                  &cumulative_rolling_resistance_spring_torque)

{
  // Useful value used more than once:
  // mu_r * R_e
  const double mu_r_times_R_e =
    effective_rolling_friction_coefficient * effective_r;


  // For calculation of rolling resistance torque, we need to obtain
  // omega_ij using rotational velocities of particles one and two
  // omega_ij is the relative angular velocity
  Tensor<1, 3> omega_ij;
  for (int d = 0; d < 3; ++d)
    // particle_one - particle_two
    omega_ij[d] = particle_one_properties[PropertiesIndex::omega_x + d] -
                  particle_two_properties[PropertiesIndex::omega_x + d];

  // Non-collinear component of the relative velocity.
  const Tensor<1, 3> omega_ij_perpendicular =
    omega_ij -
    scalar_product(omega_ij, normal_unit_vector) * normal_unit_vector;

  // Delta theta : incremental relative rotation between i and j
  const Tensor<1, 3> delta_theta = dt * omega_ij_perpendicular;

  // Rolling stiffness
  const double K_r = [&]() {
    if constexpr (dim == 3)
      return 2.25 * normal_spring_constant *
             Utilities::fixed_power<2>(mu_r_times_R_e);
    else
      return 3. * normal_spring_constant *
             Utilities::fixed_power<2>(mu_r_times_R_e);
  }();

  // Update the spring torque
  cumulative_rolling_resistance_spring_torque -= K_r * delta_theta;

  // Limiting spring torque
  const double M_r_max = mu_r_times_R_e * normal_force_norm;

  const double rolling_resistance_spring_torque_norm =
    cumulative_rolling_resistance_spring_torque.norm();

  // 1.4 * m_i * (R_i)^2
  // Total inertia of particle i evaluated at the contact point.
  const double I_i = 1.4 * particle_one_properties[PropertiesIndex::mass] *
                     Utilities::fixed_power<2>(
                       0.5 * particle_one_properties[PropertiesIndex::dp]);

  // 1.4 * m_j * (R_j)^2
  // Total inertia of particle evaluated at the contact point.
  const double I_j = 1.4 * particle_two_properties[PropertiesIndex::mass] *
                     Utilities::fixed_power<2>(
                       0.5 * particle_two_properties[PropertiesIndex::dp]);

  // Harmonic mean of the inertia at the point of contact. (Effective
  // inertia)
  const double I_e = I_i * I_j / (I_i + I_j);

  // C_r_crit = 2. * sqrt(I_r * K_r)
  // C_r = eta_r * C_r_crit
  const double C_r =
    effective_rolling_viscous_damping_coefficient * 2. * sqrt(I_e * K_r);

  // Similarly to the coulomb limit, the spring torque must be decrease to the
  // limit value if it exceeds the limiting spring toque.
  if (rolling_resistance_spring_torque_norm > M_r_max)
    {
      cumulative_rolling_resistance_spring_torque =
        cumulative_rolling_resistance_spring_torque *
        (M_r_max / rolling_resistance_spring_torque_norm);

      // If the limiting spring torque is exceeded, the damping is multiplied by
      // the f_coefficient which should be between 0 and 1.

      return cumulative_rolling_resistance_spring_torque -
             f_coefficient * C_r * omega_ij_perpendicular;
    }
  // Minus sign has been verified
  return cumulative_rolling_resistance_spring_torque -
         C_r * omega_ij_perpendicular;
}
#endif
