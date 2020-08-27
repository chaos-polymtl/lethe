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

#ifndef PWNONLINEARFORCE_H_
#  define PWNONLINEARFORCE_H_

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
  PWNonLinearForce()
  {}

  /**
   * Carries out the calculation of the particle-wall contact force using
   * non-linear (Hertzian) model
   *
   * @param pw_pairs_in_contact Required information for calculation of the
   * particle-wall contact force, these information were obtained in the
   * fine search
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  calculate_pw_contact_force(
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      *                             pw_pairs_in_contact,
    const DEMSolverParameters<dim> &dem_parameters,
    const double &                  dt) override;

private:
  /**
   * Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param physical_properties Physical properties of the system
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param particle_properties Properties of particle one in contact
   */
  std::tuple<Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>, Tensor<1, dim>>
  calculate_nonlinear_contact_force_and_torque(
    const Parameters::Lagrangian::PhysicalProperties &physical_properties,
    pw_contact_info_struct<dim> &                     contact_info,
    const ArrayView<const double> &                   particle_properties);
};

#endif
