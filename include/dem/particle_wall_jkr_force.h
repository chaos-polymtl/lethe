// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_jkr_force_h
#define lethe_particle_wall_jkr_force_h

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_contact_force.h>

#include <deal.II/particles/particle.h>

#include <cmath>
#include <iostream>
#include <vector>

using namespace dealii;

/**
 * @brief Calculation of the JKR particle-wall contact force using the
 * information obtained from the fine search and physical properties of
 * particles and walls
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticleWallJKRForce
  : public ParticleWallContactForce<dim, PropertiesIndex>
{
  using FuncPtrType =
    Tensor<1, 3> (ParticleWallJKRForce<dim, PropertiesIndex>::*)(
      const ArrayView<const double> &,
      const double,
      const double,
      const double,
      const double,
      const double,
      const Tensor<1, 3> &,
      Tensor<1, 3> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  ParticleWallJKRForce(
    const DEMSolverParameters<dim>        &dem_parameters,
    const std::vector<types::boundary_id> &boundary_index = {});

  /**
   * @brief Carries out the calculation of the particle-wall contact force using
   * JKR model
   *
   * @param[in] particle_wall_pairs_in_contact Required information for the
   * calculation of the particle-wall contact force. These information were
   * obtained in the fine search
   * @param[in] dt DEM time step
   * @param[out] torque Torque acting on particles
   * @param[out] force Force acting on particles
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                              &particle_wall_pairs_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force) override;

  /**
   * @brief Carries out the calculation of particle-floating mesh contact force using
   * JKR model
   *
   * @param[in] particle_floating_mesh_in_contact A container that stores the
   * information of particle-floating mesh contact
   * @param[in] dt DEM time step
   * @param[out] torque Torque acting on particles
   * @param[out] force Force acting on particles
   * @param[in] solids Floating solids
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                              &particle_floating_mesh_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids)
    override;


private:
  /**
   * @brief No rolling resistance torque model.
   *
   * @param[in] particle_properties Particle one properties.
   * @param[in] effective_rolling_friction_coefficient Effective rolling
   * friction coefficient.
   * @param[in] effective_rolling_viscous_damping_coefficient Effective rolling
   * viscous damping coefficient.
   * @param[in] normal_force_norm Normal force norm.
   * @param[in] dt DEM time step.
   * @param[in] normal_spring_constant Normal spring constant.
   * @param[in] normal_unit_vector Normal unit contact vector between the
   * particle and the wall.
   * @param[in,out] cumulative_rolling_resistance_spring_torque cumulative
   * rolling resistance torque for the EPSD rolling resistance model.
   *
   * @return Rolling resistance torque equal to the null vector.
   */
  inline Tensor<1, 3>
  no_resistance(const ArrayView<const double> & /*particle_properties*/,
                const double /*effective_rolling_friction_coefficient*/,
                const double /*effective_rolling_viscous_damping_coefficient*/,
                const double /*normal_force_norm*/,
                const double /*dt*/,
                const double /*normal_spring_constant*/,
                const Tensor<1, 3> & /*normal_unit_vector*/,
                Tensor<1, 3> & /*cumulative_rolling_resistance_spring_torque*/)
  {
    Tensor<1, 3> rolling_resistance({0, 0, 0});
    return rolling_resistance;
  }

  /**
   * @brief Calculation of the rolling resistance torque model using the
   * physical properties and the constant model.
   *
   * @param[in] particle_properties Particle one properties.
   * @param[in] effective_rolling_friction_coefficient Effective rolling
   * friction coefficient.
   * @param[in] effective_rolling_viscous_damping_coefficient Effective rolling
   * viscous damping coefficient.
   * @param[in] normal_force_norm Normal force norm.
   * @param[in] dt DEM time step.
   * @param[in] normal_spring_constant Normal spring constant.
   * @param[in] normal_unit_vector Normal unit contact vector between the
   * particle and the wall.
   * @param[in,out] cumulative_rolling_resistance_spring_torque cumulative
   * rolling resistance torque for the EPSD rolling resistance model.
   *
   * @return Rolling resistance torque
   */
  inline Tensor<1, 3>
  constant_resistance(
    const ArrayView<const double> &particle_properties,
    const double                   effective_rolling_friction_coefficient,
    const double /*effective_rolling_viscous_damping_coefficient */,
    const double normal_force_norm,
    const double /*dt*/,
    const double /*normal_spring_constant*/,
    const Tensor<1, 3> & /*normal_unit_vector*/,
    Tensor<1, 3> & /*cumulative_rolling_resistance_spring_torque*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] = particle_properties[PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      (particle_properties[PropertiesIndex::dp] * 0.5) * normal_force_norm *
      particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Calculation of the rolling resistance torque model using the
   * physical properties and the viscous model.
   *
   * @param particle_properties Particle one properties.
   * @param effective_rolling_friction_coefficient Effective rolling friction
   * coefficient.
   * @param effective_rolling_viscous_damping_coefficient Effective rolling
   * viscous damping coefficient.
   * @param normal_force_norm Normal force norm.
   * @param dt DEM time step.
   * @param normal_spring_constant Normal spring constant.
   * @param normal_unit_vector Normal unit contact vector between the particle
   * and the wall.
   * @param cumulative_rolling_resistance_spring_torque cumulative rolling
   * resistance torque for the EPSD rolling resistance model.
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  viscous_resistance(
    const ArrayView<const double> &particle_properties,
    const double                   effective_rolling_friction_coefficient,
    const double /*effective_rolling_viscous_damping_coefficient*/,
    const double normal_force_norm,
    const double /*dt*/,
    const double /*normal_spring_constant*/,
    const Tensor<1, 3> &normal_unit_vector,
    Tensor<1, 3> & /*cumulative_rolling_resistance_spring_torque*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] = particle_properties[PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0.)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    Tensor<1, 3> v_omega =
      cross_product_3d(angular_velocity,
                       particle_properties[PropertiesIndex::dp] * 0.5 *
                         normal_unit_vector);

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      particle_properties[PropertiesIndex::dp] * 0.5 * normal_force_norm *
      v_omega.norm() * particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief alculation of the rolling resistance torque model using the
   * physical properties and the elastic-plastic spring-dashpot model
   *
   * @param particle_properties Particle one properties.
   * @param effective_rolling_friction_coefficient Effective rolling friction
   * coefficient.
   * @param effective_rolling_viscous_damping_coefficient Effective rolling viscous
   * damping coefficient.
   * @param normal_force_norm Normal force norm.
   * @param dt DEM time step.
   * @param normal_spring_constant Normal spring constant.
   * @param normal_unit_vector Normal unit contact vector between the particle
   * and the wall.
   * @param cumulative_rolling_resistance_spring_torque cumulative rolling
   * resistance torque for the EPSD rolling resistance model.
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  epsd_resistance(const ArrayView<const double> &particle_properties,
                  const double effective_rolling_friction_coefficient,
                  const double effective_rolling_viscous_damping_coefficient,
                  const double normal_force_norm,
                  const double dt,
                  const double normal_spring_constant,
                  const Tensor<1, 3> &normal_unit_vector,
                  Tensor<1, 3> &cumulative_rolling_resistance_spring_torque)
  {
    // Useful value used more than once
    const double mu_r_times_R_e = effective_rolling_friction_coefficient * 0.5 *
                                  particle_properties[PropertiesIndex::dp];

    // Getting the angular velocity of particle in the vector format.
    Tensor<1, 3>
      omega_ij; // We ignore the rotation of the wall. j is the particle.
    for (int d = 0; d < 3; ++d)
      {
        omega_ij[d] = particle_properties[PropertiesIndex::omega_x + d];
      }

    // Non-collinear component of the relative velocity.
    const Tensor<1, 3> omega_ij_perpendicular =
      omega_ij -
      scalar_product(omega_ij, normal_unit_vector) * normal_unit_vector;

    // Delta theta
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

    // Total inertia of particle i evaluated at its surface. Mass of the
    // wall is considered infinite, thus I_i is equal to I_e. (Effective
    // inertia)
    const double I_e =
      1.4 * particle_properties[PropertiesIndex::mass] *
      Utilities::fixed_power<2>(0.5 * particle_properties[PropertiesIndex::dp]);

    // C_r_crit = 2. * sqrt(I_r * K_r)
    const double C_r =
      effective_rolling_viscous_damping_coefficient * 2. * sqrt(I_e * K_r);

    // Similarly to the coulomb limit, the spring torque must be decrease to the
    // limit value if it exceeds the limiting spring toque.
    if (rolling_resistance_spring_torque_norm > M_r_max)
      {
        cumulative_rolling_resistance_spring_torque =
          cumulative_rolling_resistance_spring_torque *
          (M_r_max / rolling_resistance_spring_torque_norm);

        // If the limiting spring torque is exceeded, the damping is multiplied
        // by the f_coefficient which should be between 0 and 1.

        return cumulative_rolling_resistance_spring_torque -
               f_coefficient_epsd * C_r * omega_ij_perpendicular;
      }
    // Minus sign has been verified
    return cumulative_rolling_resistance_spring_torque -
           C_r * omega_ij_perpendicular;
  }

  /**
   * @brief Carries out the calculation of the particle-wall JKR contact
   * force and torques based on the updated values in contact_info
   *
   * @param[in] dt DEM time step.
   * @param[in,out] contact_info A container that contains the required
   * information for calculation of the contact force for a particle pair in
   * contact.
   * @param[in] particle_properties Properties of particle one in contact
   *
   * @return A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair.
   */
  std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
  calculate_jkr_contact_force_and_torque(
    const double                     dt,
    particle_wall_contact_info<dim> &contact_info,
    const ArrayView<const double>   &particle_properties);

  // Model parameter
  const double f_coefficient_epsd;
};

#endif
