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

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_info_struct.h>

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#  define PWCONTACTFORCE_H_

/**
 * Base interface for classes that carry out the calculation of particle-wall
 * contact force
 */

template <int dim>
class PWContactForce
{
public:
  PWContactForce()
  {}

  virtual ~PWContactForce()
  {}

  /**
   * Carries out the calculation of the particle-wall contact force using the
   * contact pair information obtained from the particle-wall fine search and
   * physical properties of particles and walls
   *
   * @param pw_pairs_in_contact Required information for calculation of the
   * particle-wall contact force
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  calculate_pw_contact_force(
    std::map<int, std::map<int, pw_contact_info_struct<dim>>>
      *                             pw_pairs_in_contact,
    const DEMSolverParameters<dim> &dem_parameters,
    const double &                  dt) = 0;

protected:
  /**
   * Carries out updating the contact pair information for both non-linear and
   * linear contact force calculations
   *
   * @param contact_pair_information Contact information of a particle-wall pair
   * in neighborhood
   * @param particle_properties Properties of particle in contact
   * @param dt DEM time step
   */
  void
  update_contact_information(
    pw_contact_info_struct<dim> &  contact_pair_information,
    const ArrayView<const double> &particle_properties,
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
  apply_force_and_torque(ArrayView<double> &               particle_properties,
                         const std::tuple<Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>> &forces_and_torques);

  /** This function is used to find the projection of vector_a on
   * vector_b
   * @param vector_a A vector which is going to be projected on vector_b
   * @param vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  Tensor<1, dim> find_projection(Tensor<1, dim> vector_a,
                                 Tensor<1, dim> vector_b);
};

#endif /* PWCONTACTFORCE_H_ */
