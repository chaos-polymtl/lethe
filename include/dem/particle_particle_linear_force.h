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

#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_contact_force.h>
#include <dem/particle_particle_contact_info_struct.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

#include <math.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_particle_linear_force_h
#  define particle_particle_linear_force_h

/**
 * Calculation of the linear particle-particle contact force using the
 * information obtained from the fine search and physical properties of
 * particles
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleParticleLinearForce : public ParticleParticleContactForce<dim>
{
public:
  ParticleParticleLinearForce<dim>(
    const DEMSolverParameters<dim> &dem_parameters);

  /**
   * Carries out the calculation of the particle-particle contact force
   * using linear (Hookean) model
   *
   * @param local_adjacent_particles Required information for the calculation
   * of the local-local particle-particle contact force. This information
   * was obtained in the fine search
   * @param ghost_adjacent_particles Required information for the calculation
   * of the local-ghost particle-particle contact force. This information
   * was obtained in the fine search
   * @param dt DEM time-step
   * @param torque An unordered_map of torque of particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_particle_particle_contact_force(
    std::unordered_map<
      types::particle_index,
      std::unordered_map<types::particle_index,
                         particle_particle_contact_info_struct<dim>>>
      &local_adjacent_particles,
    std::unordered_map<
      types::particle_index,
      std::unordered_map<types::particle_index,
                         particle_particle_contact_info_struct<dim>>>
      &                        ghost_adjacent_particles,
    const double &             dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force) override;

  /**
   * Carries out the calculation of the contact force for IB particles. This
   * function is used in fem-dem/ib_particles_dem.
   *
   * @param normal_overlap Contact normal overlap. This is already calculated and
   * will be used here to calculate contact force.
   * @param contact_info Contact history including tangential overlap and relative
   * velocity.
   * @param normal_force Contact normal force.
   * @param tangential_force Contact tangential force.
   * @param particle_one_tangential_torque
   * @param particle_two_tangential_torque
   * @param rolling_resistance_torque Contact rolling resistance torque.
   * @param particle_one
   * @param particle_two
   * @param particle_one_location Location of particle one.
   * @param particle_two_location Location of particle two.
   * @param dt Time-step.
   * @param particle_one_radius radius of particle one.
   * @param particle_two_radius radius of particle two.
   * @param particle_one_mass mass of particle two.
   * @param particle_two_mass mass of particle two.
   */
  virtual void
  calculate_IB_particle_particle_contact_force(
    const double &                              normal_overlap,
    particle_particle_contact_info_struct<dim> &contact_info,
    Tensor<1, 3> &                              normal_force,
    Tensor<1, 3> &                              tangential_force,
    Tensor<1, 3> &                              particle_one_tangential_torque,
    Tensor<1, 3> &                              particle_two_tangential_torque,
    Tensor<1, 3> &                              rolling_resistance_torque,
    IBParticle<dim> &                           particle_one,
    IBParticle<dim> &                           particle_two,
    const Point<dim> &                          particle_one_location,
    const Point<dim> &                          particle_two_location,
    const double &                              dt,
    const double &                              particle_one_radius,
    const double &                              particle_two_radius,
    const double &                              particle_one_mass,
    const double &                              particle_two_mass) override;

private:
  /**
   * Carries out the calculation of the particle-particle linear contact
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
  calculate_linear_contact_force_and_torque(
    particle_particle_contact_info_struct<dim> &contact_info,
    const double &                              normal_relative_velocity_value,
    const Tensor<1, 3> &                        normal_unit_vector,
    const double &                              normal_overlap,
    const ArrayView<const double> &             particle_one_properties,
    const ArrayView<const double> &             particle_two_properties,
    Tensor<1, 3> &                              normal_force,
    Tensor<1, 3> &                              tangential_force,
    Tensor<1, 3> &                              particle_one_tangential_torque,
    Tensor<1, 3> &                              particle_two_tangential_torque,
    Tensor<1, 3> &                              rolling_resistance_torque);

  // Normal and tangential contact forces, tangential and rolling torques,
  // normal unit vector of the contact and contact relative velocity in the
  // normal direction
  Tensor<1, 3>                 normal_unit_vector;
  Tensor<1, 3>                 normal_force;
  Tensor<1, 3>                 tangential_force;
  Tensor<1, 3>                 particle_one_tangential_torque;
  Tensor<1, 3>                 particle_two_tangential_torque;
  Tensor<1, 3>                 rolling_resistance_torque;
  double                       normal_relative_velocity_value;
  RollingResistanceTorqueModel rolling_reistance_model;
};

#endif
