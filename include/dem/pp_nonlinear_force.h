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

#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_force.h>
#include <dem/pp_contact_info_struct.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

#include <math.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_particle_nonlinear_force_h
#  define particle_particle_nonlinear_force_h

/**
 * Calculation of the non-linear particle-particle contact force using the
 * information obtained from the fine search and physical properties of
 * particles
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class PPNonLinearForce : public PPContactForce<dim>
{
  using FuncPtrType =
    Tensor<1, dim> (PPNonLinearForce<dim>::*)(const double &,
                                              const ArrayView<const double> &,
                                              const ArrayView<const double> &,
                                              const double &,
                                              const double &,
                                              const Tensor<1, dim> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  PPNonLinearForce<dim>(const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Carries out the calculation of the particle-particle contact force using
   * non-linear (Hertzian) model
   *
   * @param adjacent_particles Required information for the calculation of the
   * local-local particle-particle contact force. These information were
   * obtained in the fine search
   * @param ghost_adjacent_particles Required information for the calculation of the
   * local-ghost particle-particle contact force. These information were
   * obtained in the fine search
   * @param dt DEM time-step
   * @param momentum An unordered_map of momentum of particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_pp_contact_force(
    std::unordered_map<
      types::particle_index,
      std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
      &adjacent_particles,
    std::unordered_map<
      types::particle_index,
      std::unordered_map<types::particle_index, pp_contact_info_struct<dim>>>
      &                          ghost_adjacent_particles,
    const double &               dt,
    std::vector<Tensor<1, dim>> &momentum,
    std::vector<Tensor<1, dim>> &force) override;

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
  no_resistance(const double & /*effective_r*/,
                const ArrayView<const double> & /*particle_one_properties*/,
                const ArrayView<const double> & /*particle_two_properties*/,
                const double & /*effective_rolling_friction_coefficient*/,
                const double & /*normal_force_norm*/,
                const Tensor<1, dim> & /*normal_contact_vector*/)
  {
    Tensor<1, dim> rolling_resistance;
    for (int d = 0; d < dim; ++d)
      rolling_resistance[d] = 0;

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
  constant_resistance(const double &                 effective_r,
                      const ArrayView<const double> &particle_one_properties,
                      const ArrayView<const double> &particle_two_properties,
                      const double &effective_rolling_friction_coefficient,
                      const double &normal_force_norm,
                      const Tensor<1, dim> & /*normal_contact_vector*/)
  {
    // For calculation of rolling resistance torque, we need to obtain
    // omega_ij using rotational velocities of particles one and two
    Tensor<1, dim> particle_one_angular_velocity, particle_two_angular_velocity,
      omega_ij, omega_ij_direction;
    for (int d = 0; d < dim; ++d)
      {
        particle_one_angular_velocity[d] =
          particle_one_properties[DEM::PropertiesIndex::omega_x + d];
        particle_two_angular_velocity[d] =
          particle_two_properties[DEM::PropertiesIndex::omega_x + d];
      }

    omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
    omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

    // Calculation of rolling resistance torque
    Tensor<1, dim> rolling_resistance_torque =
      -effective_rolling_friction_coefficient * effective_r *
      normal_force_norm * omega_ij_direction;

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
  viscous_resistance(const double &                 effective_r,
                     const ArrayView<const double> &particle_one_properties,
                     const ArrayView<const double> &particle_two_properties,
                     const double &effective_rolling_friction_coefficient,
                     const double &normal_force_norm,
                     const Tensor<1, dim> &normal_contact_vector)
  {
    // For calculation of rolling resistance torque, we need to obtain
    // omega_ij using rotational velocities of particles one and two
    Tensor<1, dim> particle_one_angular_velocity, particle_two_angular_velocity,
      omega_ij, omega_ij_direction;
    for (int d = 0; d < dim; ++d)
      {
        particle_one_angular_velocity[d] =
          particle_one_properties[DEM::PropertiesIndex::omega_x + d];
        particle_two_angular_velocity[d] =
          particle_two_properties[DEM::PropertiesIndex::omega_x + d];
      }

    omega_ij = particle_one_angular_velocity - particle_two_angular_velocity;
    omega_ij_direction = omega_ij / (omega_ij.norm() + DBL_MIN);

    Tensor<1, dim> v_omega =
      cross_product_3d(particle_one_angular_velocity,
                       particle_one_properties[DEM::PropertiesIndex::dp] * 0.5 *
                         normal_contact_vector) -
      cross_product_3d(particle_two_angular_velocity,
                       particle_two_properties[DEM::PropertiesIndex::dp] * 0.5 *
                         -normal_contact_vector);

    // Calculation of rolling resistance torque
    Tensor<1, dim> rolling_resistance_torque =
      -effective_rolling_friction_coefficient * effective_r *
      normal_force_norm * v_omega.norm() * omega_ij_direction;

    return rolling_resistance_torque;
  }

  /**
   * @brief Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param normal_relative_velocity_value Normal relative contact velocity
   * @param normal_unit_vector Contact normal unit vector
   * @param normal_overlap Contact normal overlap
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param particle_one_tangential_torque Contact tangential torque on particle one
   * @param particle_two_tangential_torque Contact tangential torque on particle two
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  void
  calculate_nonlinear_contact_force_and_torque(
    pp_contact_info_struct<dim> &  contact_info,
    const double &                 normal_relative_velocity_value,
    const Tensor<1, dim> &         normal_unit_vector,
    const double &                 normal_overlap,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_propertie,
    Tensor<1, dim> &               normal_force,
    Tensor<1, dim> &               tangential_force,
    Tensor<1, dim> &               particle_one_tangential_torque,
    Tensor<1, dim> &               particle_two_tangential_torque,
    Tensor<1, dim> &               rolling_resistance_torque);

  // Contact model parameter. It is calculated in the constructor for different
  // combinations of particle types. For different combinations, a map of map is
  // used to store this variable
  std::map<int, std::map<int, double>> model_parameter_beta;

  // Normal and tangential contact forces, tangential and rolling torques and
  // normal unit vector of the contact
  Tensor<1, dim> normal_unit_vector;
  Tensor<1, dim> normal_force;
  Tensor<1, dim> tangential_force;
  Tensor<1, dim> particle_one_tangential_torque;
  Tensor<1, dim> particle_two_tangential_torque;
  Tensor<1, dim> rolling_resistance_torque;
  double         normal_relative_velocity_value;
};

#endif
