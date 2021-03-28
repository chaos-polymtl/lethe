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

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */

#include <deal.II/particles/particle.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_force.h>
#include <dem/pw_contact_info_struct.h>
#include <math.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_wall_linear_force_h
#  define particle_wall_linear_force_h

/**
 * Calculation of the linear particle-wall contact force using the
 * information obtained from the fine search and physical properties of
 * particles and walls
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class PWLinearForce : public PWContactForce<dim>
{
  using FuncPtrType =
    Tensor<1, dim> (PWLinearForce<dim>::*)(const ArrayView<const double> &,
                                           const double &,
                                           const double &,
                                           const Tensor<1, dim> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  PWLinearForce<dim>(
    const std::unordered_map<types::particle_index, Tensor<1, dim>>
      boundary_translational_velocity,
    const std::unordered_map<types::particle_index, double>
      boundary_rotational_speed,
    const std::unordered_map<types::particle_index, Tensor<1, dim>>
                                    boundary_rotational_vector,
    const double                    triangulation_radius,
    const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Carries out the calculation of the particle-wall contact force using
   * linear (Hookean) model
   *
   * @param pw_pairs_in_contact Required information for calculation of
   * the particle-wall contact force. These information were obtained in
   * the fine search
   * @param dt DEM time step
   * @param momentum An unordered_map of momentum of particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_pw_contact_force(
    std::unordered_map<
      types::particle_index,
      std::map<types::particle_index, pw_contact_info_struct<dim>>>
      &           pw_pairs_in_contact,
    const double &dt,
    std::unordered_map<types::particle_index, Tensor<1, dim>> &momentum,
    std::unordered_map<types::particle_index, Tensor<1, dim>> &force) override;

private:
  /**
   * @brief No rolling resistance torque model
   *
   * @param particle_one_properties Particle one properties
   * @param particle_two_properties Particle two properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, dim>
  no_resistance(const ArrayView<const double> & /*particle_properties*/,
                const double & /*effective_rolling_friction_coefficient*/,
                const double & /*normal_force_norm*/,
                const Tensor<1, dim> & /*normal_contact_vector*/)
  {
    Tensor<1, dim> rolling_resistance;
    rolling_resistance[0] = 0;
    rolling_resistance[1] = 0;
    rolling_resistance[2] = 0;

    return rolling_resistance;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the constant model
   *
   * @param particle_one_properties Particle one properties
   * @param particle_two_properties Particle two properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, dim>
  constant_resistance(const ArrayView<const double> &particle_properties,
                      const double &effective_rolling_friction_coefficient,
                      const double &normal_force_norm,
                      const Tensor<1, dim> & /*normal_contact_vector*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, dim> angular_velocity;
    for (int d = 0; d < dim; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, dim> pw_angular_velocity;
    for (int d = 0; d < dim; ++d)
      {
        pw_angular_velocity[d] = 0;
      }

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        pw_angular_velocity = angular_velocity / omega_value;
      }

    // Calcualation of rolling resistance torque
    Tensor<1, dim> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      (particle_properties[DEM::PropertiesIndex::dp] * 0.5) *
      normal_force_norm * pw_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the viscous model
   *
   * @param particle_one_properties Particle one properties
   * @param particle_two_properties Particle two properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, dim>
  viscous_resistance(const ArrayView<const double> &particle_properties,
                     const double &effective_rolling_friction_coefficient,
                     const double &normal_force_norm,
                     const Tensor<1, dim> &normal_contact_vector)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, dim> angular_velocity;
    for (int d = 0; d < dim; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, dim> pw_angular_velocity;
    for (int d = 0; d < dim; ++d)
      {
        pw_angular_velocity[d] = 0;
      }

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        pw_angular_velocity = angular_velocity / omega_value;
      }

    Tensor<1, dim> v_omega =
      cross_product_3d(angular_velocity,
                       particle_properties[DEM::PropertiesIndex::dp] * 0.5 *
                         normal_contact_vector);

    // Calculation of rolling resistance torque
    Tensor<1, dim> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      particle_properties[DEM::PropertiesIndex::dp] * 0.5 * normal_force_norm *
      v_omega.norm() * pw_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * Carries out the calculation of the particle-particle linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param particle_properties Properties of particle one in contact
   * @return A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair
   */
  std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
  calculate_linear_contact_force_and_torque(
    pw_contact_info_struct<dim> &  contact_info,
    const ArrayView<const double> &particle_properties);
};

#endif
