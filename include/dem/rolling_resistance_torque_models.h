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
 * Author: Bruno Blais, Polytechnique Montreal, 2019
 */

#include <dem/dem_solver_parameters.h>
#include <dem/rolling_resistance_torque_base.h>


#ifndef rolling_resistance_torque_models_h
#  define rolling_resistance_torque_models_h

template <int dim>
class NoRollingResistanceTorque : public RollingResistanceTorqueBase<dim>
{
public:
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
  virtual Tensor<1, dim>
  calculate_rolling_resistance_torque(
    const double &                 effective_r,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const double &                 effective_rolling_friction_coefficient,
    const double &                 normal_force_norm,
    const Tensor<1, dim> &         normal_contact_vector) override;
};

template <int dim>
class ConstantRollingResistanceTorque : public RollingResistanceTorqueBase<dim>
{
public:
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
  virtual Tensor<1, dim>
  calculate_rolling_resistance_torque(
    const double &                 effective_r,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const double &                 effective_rolling_friction_coefficient,
    const double &                 normal_force_norm,
    const Tensor<1, dim> &         normal_contact_vector) override;
};

template <int dim>
class ViscousRollingResistanceTorque : public RollingResistanceTorqueBase<dim>
{
public:
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
  virtual Tensor<1, dim>
  calculate_rolling_resistance_torque(
    const double &                 effective_r,
    const ArrayView<const double> &particle_one_properties,
    const ArrayView<const double> &particle_two_properties,
    const double &                 effective_rolling_friction_coefficient,
    const double &                 normal_force_norm,
    const Tensor<1, dim> &         normal_contact_vector) override;
};

#endif
