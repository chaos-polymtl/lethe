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

#include <deal.II/particles/particle_handler.h>

#include <boost/range/adaptor/map.hpp>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_info_struct.h>

using namespace dealii;

#ifndef particle_particle_contact_force_h
#  define particle_particle_contact_force_h

/**
 * Base interface for classes that carry out the calculation of particle-paricle
 * contact force
 */
template <int dim>
class PPContactForce
{
public:
  PPContactForce()
  {}

  virtual ~PPContactForce()
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
   */
  virtual void
  calculate_pp_contact_force(
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
      &local_adjacent_particles,
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
      &           ghost_adjacent_particles,
    const double &dt) = 0;

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
    pp_contact_info_struct<dim> &  adjacent_pair_information,
    double &                       normal_relative_velocity_value,
    Tensor<1, dim> &               normal_unit_vector,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const Point<dim> &             particle_one_location,
    const Point<dim> &             particle_two_location,
    const double &                 dt);

  /**
   * @brief Carries out applying the calculated force and torque on the local-local
   * particle pair in contact, for both non-linear and linear contact force
   * calculations
   *
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  void
  apply_force_and_torque_real(ArrayView<double> &   particle_one_properties,
                              ArrayView<double> &   particle_two_properties,
                              const Tensor<1, dim> &normal_force,
                              const Tensor<1, dim> &tangential_force,
                              const Tensor<1, dim> &tangential_torque,
                              const Tensor<1, dim> &rolling_resistance_torque);

  /**
   * Carries out applying the calculated force and torque on the local-ghost
   * particle pair in contact, for both non-linear and linear contact force
   * calculations. The contact force is only applied on the local particles
   *
   * @param particle_one_properties Properties of particle one (local) in
   * contact
   * @param normal_force normal_force Contact normal force
   * @param tangential_force Contact tangential force
   * @param tangential_torque Contact tangential torque
   * @param rolling_friction_torque Contact rolling resistance torque
   */
  void
  apply_force_and_torque_ghost(ArrayView<double> &   particle_one_properties,
                               const Tensor<1, dim> &normal_force,
                               const Tensor<1, dim> &tangential_force,
                               const Tensor<1, dim> &tangential_torque,
                               const Tensor<1, dim> &rolling_resistance_torque);

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
