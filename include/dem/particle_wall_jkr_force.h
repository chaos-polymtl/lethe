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
 */
template <int dim, DEM::SolverType solver_type>
class ParticleWallJKRForce : public ParticleWallContactForce<dim, solver_type>
{
  using FuncPtrType = Tensor<1, 3> (ParticleWallJKRForce<dim, solver_type>::*)(
    const ArrayView<const double> &,
    const double,
    const double,
    const Tensor<1, 3> &);
  FuncPtrType calculate_rolling_resistance_torque;

public:
  ParticleWallJKRForce(
    const DEMSolverParameters<dim>        &dem_parameters,
    const std::vector<types::boundary_id> &boundary_index = {});

  /**
   * @brief Carries out the calculation of the particle-wall contact force using
   * JKR model
   *
   * @param particle_wall_pairs_in_contact Required information for the calculation of
   * the particle-wall contact force. These information were obtained in
   * the fine search
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
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
   * @param particle_floating_mesh_in_contact A container that stores the information of
   * particle-floating mesh contact
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param solids Floating solids
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
   * @brief No rolling resistance torque model
   *
   * @param particle_one_properties Particle one properties
   * @param particle_two_properties Particle two properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return Rolling resistance torque
   */
  inline Tensor<1, 3>
  no_resistance(const ArrayView<const double> & /*particle_properties*/,
                const double /*effective_rolling_friction_coefficient*/,
                const double /*normal_force_norm*/,
                const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    Tensor<1, 3> rolling_resistance({0, 0, 0});
    return rolling_resistance;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the constant model
   *
   * @param particle_properties Particle properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return Rolling resistance torque
   */
  inline Tensor<1, 3>
  constant_resistance(const ArrayView<const double> &particle_properties,
                      const double effective_rolling_friction_coefficient,
                      const double normal_force_norm,
                      const Tensor<1, 3> & /*normal_contact_vector*/)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
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
      (particle_properties[DEM::PropertiesIndex<solver_type>::dp] * 0.5) *
      normal_force_norm * particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Carries out calculation of the rolling resistance torque using the viscous model
   *
   * @param particle_properties Particle properties
   * @param effective_rolling_friction_coefficient Effective rolling friction coefficient
   * @param normal_force_norm Normal force norm
   *
   * @return rolling resistance torque
   */
  inline Tensor<1, 3>
  viscous_resistance(const ArrayView<const double> &particle_properties,
                     const double        effective_rolling_friction_coefficient,
                     const double        normal_force_norm,
                     const Tensor<1, 3> &normal_contact_vector)
  {
    // Getting the angular velocity of particle in the vector format
    Tensor<1, 3> angular_velocity;
    for (int d = 0; d < 3; ++d)
      {
        angular_velocity[d] =
          particle_properties[DEM::PropertiesIndex<solver_type>::omega_x + d];
      }

    // Calculation of particle-wall angular velocity (norm of the
    // particle angular velocity)
    Tensor<1, 3> particle_wall_angular_velocity({0.0, 0.0, 0.0});

    double omega_value = angular_velocity.norm();
    if (omega_value != 0)
      {
        particle_wall_angular_velocity = angular_velocity / omega_value;
      }

    Tensor<1, 3> v_omega = cross_product_3d(
      angular_velocity,
      particle_properties[DEM::PropertiesIndex<solver_type>::dp] * 0.5 *
        normal_contact_vector);

    // Calculation of rolling resistance torque
    Tensor<1, 3> rolling_resistance_torque =
      -effective_rolling_friction_coefficient *
      particle_properties[DEM::PropertiesIndex<solver_type>::dp] * 0.5 *
      normal_force_norm * v_omega.norm() * particle_wall_angular_velocity;

    return rolling_resistance_torque;
  }

  /**
   * @brief Carries out the calculation of the particle-wall JKR contact
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
  calculate_jkr_contact_force_and_torque(
    particle_wall_contact_info<dim> &contact_info,
    const ArrayView<const double>   &particle_properties);
};

#endif
