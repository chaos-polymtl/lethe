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

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_contact_force.h>
#include <dem/particle_wall_contact_info_struct.h>

#include <deal.II/particles/particle.h>

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
class ParticleWallLinearForce : public ParticleWallContactForce<dim>
{
  using FuncPtrType = Tensor<1, 3> (ParticleWallLinearForce<dim>::*)(
    const ArrayView<const double> &,
    const double &,
    const double &,
    const Tensor<1, 3> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  ParticleWallLinearForce<dim>(
    const std::unordered_map<unsigned int, Tensor<1, 3>>
                                                   boundary_translational_velocity,
    const std::unordered_map<unsigned int, double> boundary_rotational_speed,
    const std::unordered_map<unsigned int, Tensor<1, 3>>
                                          boundary_rotational_vector,
    const double                          triangulation_radius,
    const DEMSolverParameters<dim> &      dem_parameters,
    const std::vector<types::boundary_id> boundary_index = {});

  /**
   * Carries out the calculation of the particle-wall contact force using
   * linear (Hookean) model
   *
   * @param particle_wall_pairs_in_contact Required information for calculation of
   * the particle-wall contact force. These information were obtained in
   * the fine search
   * @param dt DEM time step
   * @param torque Torque acting on particle
   * @param force Force acting on particles
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_in_contact &particle_wall_pairs_in_contact,
    const double                      dt,
    std::vector<Tensor<1, 3>> &       torque,
    std::vector<Tensor<1, 3>> &       force) override;

  /**
   * Carries out the calculation of particle-floating mesh contact force using
   * linear (Hooekan) model
   *
   * @param particle_floating_mesh_in_contact A container that stores the information of
   * particle-floating mesh contact
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param solids Floating solids
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename dem_data_containers::dem_data_structures<dim>::
      particle_floating_mesh_in_contact &particle_floating_mesh_in_contact,
    const double                         dt,
    std::vector<Tensor<1, 3>> &          torque,
    std::vector<Tensor<1, 3>> &          force,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids)
    override;

  /**
   * Carries out the calculation of the contact force for IB particles. This
   * function is used in fem-dem/ib_particles_dem.
   *
   * TODO : At this point in time this function is not implemented and the
   * function used for particle wall contact force is hardcoded. This will be
   * modified in futur update
   *
   * @param contact_info Contact history including tangential overlap and relative
   * velocity.
   * @param normal_force Contact normal force.
   * @param tangential_force Contact tangential force.
   * @param tangential_torque
   * @param rolling_resistance_torque Contact rolling resistance torque.
   * @param particle
   * @param wall_youngs_modulus
   * @param wall_poisson_ratio
   * @param wall_restitution_coefficient
   * @param wall_friction_coefficient
   * @param wall_rolling_friction_coefficient
   * @param dt Time-step.
   */
  virtual void
  calculate_IB_particle_wall_contact_force(
    particle_wall_contact_info_struct<dim> &contact_info,
    Tensor<1, 3> &                          normal_force,
    Tensor<1, 3> &                          tangential_force,
    Tensor<1, 3> &                          tangential_torque,
    Tensor<1, 3> &                          rolling_resistance_torque,
    IBParticle<dim> &                       particle,
    const double &                          wall_youngs_modulus,
    const double &                          wall_poisson_ratio,
    const double &                          wall_restitution_coefficient,
    const double &                          wall_friction_coefficient,
    const double &                          wall_rolling_friction_coefficient,
    const double                            dt,
    const double &                          mass,
    const double &                          radius) override;

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
  inline Tensor<1, 3>
  no_resistance(const ArrayView<const double> & /*particle_properties*/,
                const double & /*effective_rolling_friction_coefficient*/,
                const double & /*normal_force_norm*/,
                const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    Tensor<1, 3> rolling_resistance({0, 0, 0});
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
  inline Tensor<1, 3>
  constant_resistance(const ArrayView<const double> &particle_properties,
                      const double &effective_rolling_friction_coefficient,
                      const double &normal_force_norm,
                      const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    // Calcualation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      (particle_properties[DEM::PropertiesIndex::dp] * 0.5) *
      normal_force_norm * particle_wall_angular_velocity;

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
  inline Tensor<1, 3>
  viscous_resistance(const ArrayView<const double> &particle_properties,
                     const double &      effective_rolling_friction_coefficient,
                     const double &      normal_force_norm,
                     const Tensor<1, 3> &normal_contact_vector)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    Tensor<1, 3> v_omega =
      cross_product_3d(angular_velocity,
                       particle_properties[DEM::PropertiesIndex::dp] * 0.5 *
                         normal_contact_vector);

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      particle_properties[DEM::PropertiesIndex::dp] * 0.5 * normal_force_norm *
      v_omega.norm() * particle_wall_angular_velocity;

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
  std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
  calculate_linear_contact_force_and_torque(
    particle_wall_contact_info_struct<dim> &contact_info,
    const ArrayView<const double> &         particle_properties);
};

#endif
