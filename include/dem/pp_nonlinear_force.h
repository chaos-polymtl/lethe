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
#include <deal.II/particles/particle_iterator.h>

#include <dem/dem_solver_parameters.h>
#include <dem/pp_contact_force.h>
#include <dem/pp_contact_info_struct.h>
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
public:
  PPNonLinearForce()
  {}

  /**
   * Carries out the calculation of the particle-particle contact force using
   * non-linear (Hertzian) model
   *
   * @param local_adjacent_particles Required information for calculation of the
   * local-local particle-particle contact force, these information were
   * obtained in the fine search
   * @param ghost_adjacent_particles Required information for calculation of the
   * local-ghost particle-particle contact force, these information were
   * obtained in the fine search
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param dt DEM time-step
   */
  virtual void
  calculate_pp_contact_force(
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
      *adjacent_particles,
    std::unordered_map<int,
                       std::unordered_map<int, pp_contact_info_struct<dim>>>
      *                                               ghost_adjacent_particles,
    const Parameters::Lagrangian::PhysicalProperties &physical_properties,
    const double &                                    dt) override;

private:
  /**
   * Carries out the calculation of the particle-particle non-linear contact
   * force and torques based on the updated values in contact_info
   *
   * @param physical_properties Physical properties of the system
   * @param contact_info A container that contains the required information for
   * calculation of the contact force for a particle pair in contact
   * @param particle_one_properties Properties of particle one in contact
   * @param particle_two_properties Properties of particle two in contact
   */
  void
  calculate_nonlinear_contact_force_and_torque(
    const Parameters::Lagrangian::PhysicalProperties &physical_properties,
    pp_contact_info_struct<dim> &                     contact_info,
    const double &                 normal_relative_velocity_value,
    const Tensor<1, dim> &         normal_unit_vector,
    const double &                 normal_overlap,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_propertie,
    Tensor<1, dim> &               normal_force,
    Tensor<1, dim> &               tangential_force,
    Tensor<1, dim> &               tangential_torque,
    Tensor<1, dim> &               rolling_resistance_torque);

  double         normal_relative_velocity_value;
  Tensor<1, dim> normal_unit_vector;
  Tensor<1, dim> normal_force;
  Tensor<1, dim> tangential_force;
  Tensor<1, dim> tangential_torque;
  Tensor<1, dim> rolling_resistance_torque;
};

#endif
