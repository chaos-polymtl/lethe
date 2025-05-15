// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_rolling_resistance_torque_h
#define lethe_particle_wall_rolling_resistance_torque_h

/**
 * @brief No rolling resistance torque. No calculation is being done with
 * this model.
 */

inline Tensor<1, 3>
no_rolling_torque()
{
  Tensor<1, 3> rolling_resistance_torque({0.0, 0.0, 0.0});
  return rolling_resistance_torque;
}

/**
 * @brief Calculation of the constant rolling resistance torque using the
 * information obtained from the fine search and physical properties.
 *
 * The model works as follows: (particle i, wall j)
 * M_r = - mu_r * R_i * |F_n| * omega_hat
 * omega_hat = omega_i / |omega_i|
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] particle_radius Particle radius.
 * @param[in] particle_one_properties Properties of particle.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the normal force.
 *
 * @return rolling_resistance_torque with constant resistance
 */
template <typename PropertiesIndex>
inline Tensor<1, 3>
constant_rolling_torque(const double                   particle_radius,
                        const ArrayView<const double> &particle_properties,
                        const double effective_rolling_friction_coefficient,
                        const double normal_force_norm)

{
  // Getting the angular velocity of the particle
  Tensor<1, 3> particle_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_angular_velocity[d] =
        particle_properties[PropertiesIndex::omega_x + d];
    }

  // Calculation of particle-wall angular velocity (norm of the
  // particle angular velocity)
  const double omega_value = particle_angular_velocity.norm();
  Tensor<1, 3> particle_wall_angular_velocity =
    particle_angular_velocity / (omega_value + DBL_MIN);

  // Calculation of rolling resistance torque
  Tensor<1, 3> rolling_resistance_torque =
    -effective_rolling_friction_coefficient * particle_radius *
    normal_force_norm * particle_wall_angular_velocity;

  return rolling_resistance_torque;
}

/**
 * @brief Calculation of the viscous rolling resistance torque using the
 * information obtained from the fine search and physical properties.
 *
 * The model works as follows: (particle i, wall j)
 * M_r = - mu_r * R_i * |F_n| * |V_omega| * omega_hat
 * omega_hat = omega_i / |omega_i|
 * V_omega = omega_i × (R_p * n_ji)
 *
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] particle_radius Particle radius.
 * @param[in] particle_properties Properties of particle.
 * @param[in] effective_rolling_friction_coefficient Effective_rolling friction
 * coefficient
 * @param[in] normal_force_norm Norm of the normal force.
 * @param[in] normal_unit_vector Normal unit vector.
 *
 * @return rolling_resistance_torque with viscous resistance
 */

template <typename PropertiesIndex>
inline Tensor<1, 3>
viscous_rolling_torque(const double                   particle_radius,
                       const ArrayView<const double> &particle_properties,
                       const double effective_rolling_friction_coefficient,
                       const double normal_force_norm,
                       const Tensor<1, 3> &normal_unit_vector)

{
  // Getting the angular velocity of particle
  Tensor<1, 3> particle_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_angular_velocity[d] =
        particle_properties[PropertiesIndex::omega_x + d];
    }

  // Calculation of particle-wall angular velocity (norm of the
  // particle angular velocity)
  const double omega_value = particle_angular_velocity.norm();
  Tensor<1, 3> particle_wall_angular_velocity =
    particle_angular_velocity / (omega_value + DBL_MIN);

  Tensor<1, 3> v_omega = cross_product_3d(particle_angular_velocity,
                                          particle_radius * normal_unit_vector);

  // Calculation of rolling resistance torque
  Tensor<1, 3> rolling_resistance_torque =
    -effective_rolling_friction_coefficient * particle_radius *
    normal_force_norm * v_omega.norm() * particle_wall_angular_velocity;

  return rolling_resistance_torque;
}

/**
 * @brief Calculation of the EPSD rolling resistance torque using the
 * information obtained from the fine search and physical properties.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the problem is solved.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param[in] particle_radius Particle radius.
 * @param[in] particle_properties Properties of particle.
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
 *
 * @return rolling_resistance_torque with EPSD resistance
 */

template <int dim, typename PropertiesIndex>
inline Tensor<1, 3>
epsd_rolling_torque(const double                   particle_radius,
                    const ArrayView<const double> &particle_properties,
                    const double effective_rolling_friction_coefficient,
                    const double effective_rolling_viscous_damping_coefficient,
                    const double f_coefficient,
                    const double normal_force_norm,
                    const double dt,
                    const double normal_spring_constant,
                    const Tensor<1, 3> &normal_unit_vector,
                    Tensor<1, 3> &cumulative_rolling_resistance_spring_torque)

{
  // Useful value used more than once:
  // mu_r * R
  const double mu_r_times_R =
    effective_rolling_friction_coefficient * particle_radius;

  // Getting the angular velocity of particle
  Tensor<1, 3> particle_angular_velocity;
  for (int d = 0; d < 3; ++d)
    {
      particle_angular_velocity[d] =
        particle_properties[PropertiesIndex::omega_x + d];
    }

  // Non-collinear component of the relative velocity.
  const Tensor<1, 3> omega_perpendicular =
    particle_angular_velocity -
    scalar_product(particle_angular_velocity, normal_unit_vector) *
      normal_unit_vector;

  // Delta theta : incremental relative rotation between i and j
  const Tensor<1, 3> delta_theta = dt * omega_perpendicular;

  // Rolling stiffness
  const double K_r = [&]() {
    if constexpr (dim == 3)
      return 2.25 * normal_spring_constant *
             Utilities::fixed_power<2>(mu_r_times_R);
    else
      return 3. * normal_spring_constant *
             Utilities::fixed_power<2>(mu_r_times_R);
  }();

  // Update the spring torque
  cumulative_rolling_resistance_spring_torque -= K_r * delta_theta;

  // Limiting spring torque
  const double M_r_max = mu_r_times_R * normal_force_norm;

  const double rolling_resistance_spring_torque_norm =
    cumulative_rolling_resistance_spring_torque.norm();

  // 1.4 * m_i * (R_i)^2
  // Total inertia of particle i evaluated at its surface. Mass of the
  // wall is considered infinite, thus I_i is equal to I_e. (Effective
  // inertia)
  const double I_e = 1.4 * particle_properties[PropertiesIndex::mass] *
                     Utilities::fixed_power<2>(particle_radius);

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
             f_coefficient * C_r * omega_perpendicular;
    }
  // Minus sign has been verified
  return cumulative_rolling_resistance_spring_torque -
         C_r * omega_perpendicular;
}
#endif
