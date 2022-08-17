/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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

#include <core/auxiliary_math_functions.h>
#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_particle_contact_info_struct.h>
#include <dem/rolling_resistance_torque_models.h>

#include <deal.II/particles/particle_handler.h>

#include <boost/range/adaptor/map.hpp>

using namespace dealii;

#ifndef particle_particle_contact_force_h
#  define particle_particle_contact_force_h

/**
 * Base interface for classes that carry out the calculation of
 * particle-particle contact force including non-linear and linear contact
 * models
 */
template <int dim>
class ParticleParticleContactForce
{
public:
  ParticleParticleContactForce<dim>()
  {}

  virtual ~ParticleParticleContactForce()
  {}

  /**
   * Carries out the calculation of the contact force using the contact pair
   * information
   * obtained in the fine search and physical properties of particles
   *
   * @param local_adjacent_particles Required information for calculation of the
   * loacl-local particle-particle contact force. These information were
   * obtained in the fine search
   * @param ghost_adjacent_particles Required information for calculation of the
   * loacl-ghost particle-particle contact force. These information were
   * obtained in the fine search
   * @param dt DEM time step
   * @param torque An unordered_map of torque of particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_particle_particle_contact_force(
    typename dem_data_containers::dem_data_structures<
      dim>::adjacent_particle_pairs &local_adjacent_particles,
    typename dem_data_containers::dem_data_structures<
      dim>::adjacent_particle_pairs &ghost_adjacent_particles,
    const double                     dt,
    std::vector<Tensor<1, 3>> &      torque,
    std::vector<Tensor<1, 3>> &      force) = 0;

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
    const double                                normal_overlap,
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
    const double                                dt,
    const double                                particle_one_radius,
    const double                                particle_two_radius,
    const double                                particle_one_mass,
    const double                                particle_two_mass) = 0;

protected:
  /**
   * @brief Carries out updating the contact pair information for both non-linear and
   * linear contact force calculations
   *
   * @param adjacent_pair_information Contact information of a particle pair in
   * neighborhood
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param particle_one_location Location of particle one in contact
   * @param particle_two_location Location of particle two in contact
   * @param dt DEM time step
   */
  void
  update_contact_information(
    particle_particle_contact_info_struct<dim> &adjacent_pair_information,
    double &                                    normal_relative_velocity_value,
    Tensor<1, 3> &                              normal_unit_vector,
    const ArrayView<const double> &             particle_one_properties,
    const ArrayView<const double> &             particle_two_properties,
    const Point<3> &                            particle_one_location,
    const Point<3> &                            particle_two_location,
    const double                                dt);

  /**
   * @brief Carries out applying the calculated force and torque on the local-local
   * particle pair in contact, for both non-linear and linear contact force
   * calculations
   *
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   * @param particle_one_torque
   * @param particle_two_torque
   * @param particle_one_force Force acting on particle one
   * @param particle_two_force Force acting on particle two
   */
  inline void
  apply_force_and_torque_on_local_particles(
    const Tensor<1, 3> &normal_force,
    const Tensor<1, 3> &tangential_force,
    const Tensor<1, 3> &particle_one_tangential_torque,
    const Tensor<1, 3> &particle_two_tangential_torque,
    const Tensor<1, 3> &rolling_resistance_torque,
    Tensor<1, 3> &      particle_one_torque,
    Tensor<1, 3> &      particle_two_torque,
    Tensor<1, 3> &      particle_one_force,
    Tensor<1, 3> &      particle_two_force)
  {
    // Calculation of total force
    Tensor<1, 3> total_force = normal_force + tangential_force;

    // Updating the force and torque of particles in the particle handler
    particle_one_force -= total_force;
    particle_two_force += total_force;
    particle_one_torque +=
      -particle_one_tangential_torque + rolling_resistance_torque;
    particle_two_torque +=
      -particle_two_tangential_torque - rolling_resistance_torque;
  }

  /**
   * Carries out applying the calculated force and torque on the local-ghost
   * particle pair in contact, for both non-linear and linear contact force
   * calculations. The contact force is only applied on the local particles
   *
   * @param normal_force normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   * @param particle_one_torque Torque acting on particle one (local)
   * @param particle_one_force Force acting on particle one
   */
  inline void
  apply_force_and_torque_on_ghost_particles(
    const Tensor<1, 3> &normal_force,
    const Tensor<1, 3> &tangential_force,
    const Tensor<1, 3> &particle_one_tangential_torque,
    const Tensor<1, 3> &rolling_resistance_torque,
    Tensor<1, 3> &      particle_one_torque,
    Tensor<1, 3> &      particle_one_force)
  {
    // Calculation of total force
    Tensor<1, 3> total_force = normal_force + tangential_force;

    // Updating the force and torque acting on particles in the particle handler
    particle_one_force -= total_force;
    particle_one_torque +=
      -particle_one_tangential_torque + rolling_resistance_torque;
  }

  /**
   * Carries out the calculation of effective mass and radius of particles i and
   * j in contact.
   *
   * @param particle_one_properties Properties of particle one in
   * contact
   * @param particle_two_properties Properties of particle two in
   * contact
   */
  void
  find_effective_radius_and_mass(
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties);


  std::map<int, std::map<int, double>> effective_youngs_modulus;
  std::map<int, std::map<int, double>> effective_shear_modulus;
  std::map<int, std::map<int, double>> effective_coefficient_of_restitution;
  std::map<int, std::map<int, double>> effective_coefficient_of_friction;
  std::map<int, std::map<int, double>>
         effective_coefficient_of_rolling_friction;
  double effective_radius;
  double effective_mass;
};

#endif /* particle_particle_contact_force_h */
