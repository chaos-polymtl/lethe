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

#ifndef particle_wall_nonlinear_force_h
#  define particle_wall_nonlinear_force_h

/**
 * Calculation of the non-linear particle-wall contact force using the
 * information obtained from the fine search and physical properties of
 * particles and walls
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class PWNonLinearForce : public PWContactForce<dim>
{
public:
  PWNonLinearForce<dim>(
    const std::unordered_map<int, Tensor<1, dim>>
                                                  boundary_translational_velocity,
    const std::unordered_map<int, double>         boundary_rotational_speed,
    const std::unordered_map<int, Tensor<1, dim>> boundary_rotational_vector,
    const double                                  triangulation_radius,
    const DEMSolverParameters<dim> &              dem_parameters);

  /**
   * Carries out the calculation of the particle-wall contact force using
   * non-linear (Hertzian) model
   *
   * @param pw_pairs_in_contact Required information for calculation of
   * the particle-wall contact force, these information were obtained in
   * the fine search
   * @param dt DEM time step
   */
  virtual void
  calculate_pw_contact_force(
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      &           pw_pairs_in_contact,
    const double &dt) override;

private:
  /**
   * Carries out the calculation of the particle-particle non-linear contact
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
  calculate_nonlinear_contact_force_and_torque(
    pw_contact_info_struct<dim> &  contact_info,
    const ArrayView<const double> &particle_properties);
};

#endif
