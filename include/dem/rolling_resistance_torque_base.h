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

//#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/particles/particle_handler.h>


using namespace dealii;

#ifndef rolling_resistance_torque_base_h
#  define rolling_resistance_torque_base_h

/**
 * Base interface for classes that carry out the calculation of rolling
 * resistance torque
 */
template <int dim>
class RollingResistanceTorqueBase
{
public:
  RollingResistanceTorqueBase<dim>()
  {}

  virtual ~RollingResistanceTorqueBase()
  {}

  /**
   * Carries out the calculation of the rolling resistance torque
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
    const Tensor<1, dim> &         normal_contact_vector) = 0;
};

#endif /* rolling_resistance_torque_base_h */
