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

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_info_struct.h>

using namespace dealii;

#ifndef PPCONTACTFORCE_H_
#  define PPCONTACTFORCE_H_

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
   * @param pairs_in_contact_info Required information for calculation of the
   * particle-particle contact force, these information were obtained in the
   * fine search
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  calculate_pp_contact_force(
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
      *local_adjacent_particles,
    std::map<int, std::map<int, pp_contact_info_struct<dim>>>
      *                             ghost_adjacent_particles,
    const DEMSolverParameters<dim> &dem_parameters,
    const double &                  dt) = 0;

protected:
  /**
   * Carries out updating the contact pair information for both non-linear and
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
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const Point<dim> &             particle_one_location,
    const Point<dim> &             particle_two_location,
    const double &                 dt);

  /**
   * Carries out applying the calculated force and torque on the local-local
   * particle pair in contact, for both non-linear and linear contact force
   * calculations
   *
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   * @param forces_and_torques A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair
   */
  void
  apply_force_and_torque_real(
    ArrayView<double> &particle_one_properties,
    ArrayView<double> &particle_two_properties,
    const std::
      tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
        &forces_and_torques);

  /**
   * Carries out applying the calculated force and torque on the local-ghost
   * particle pair in contact, for both non-linear and linear contact force
   * calculations. The contact force is only applied on the local particles
   *
   * @param particle_one_properties Properties of particle one (local) in
   * contact
   * @param forces_and_torques A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair
   */
  void
  apply_force_and_torque_ghost(
    ArrayView<double> &particle_one_properties,
    const std::
      tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
        &forces_and_torques);
};

#endif /* PPCONTACTFORCE_H_ */
