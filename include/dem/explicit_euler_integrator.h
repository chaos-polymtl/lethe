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

#include <deal.II/particles/particle_handler.h>

#include <dem/dem_solver_parameters.h>
#include <dem/integrator.h>

using namespace dealii;

#ifndef explicit_Euler_integrator_h
#  define explicit_Euler_integrator_h

/**
 * Implementation of a classical explicit euler scheme for the integration
 * of the particle motion. Note that reinitilization of force and torque is also
 * integrated into integration class
 *
 * @note Euler is a first-order integration scheme. Calculation proceudre:
 *
 * x(n+1) = x(n) + v(n) * dt
 * v(n+1) = v(n) + a(n) * dt
 * a(n+1) = F(n+1) / m
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class ExplicitEulerIntegrator : public Integrator<dim>
{
public:
  ExplicitEulerIntegrator()
  {}

  /**
   * Carries out the prediction (pre_force) integration calculations.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   */
  virtual void
  integrate_pre_force(Particles::ParticleHandler<dim> &particle_handler,
                      Tensor<1, dim>                   body_force,
                      double                           time_step) override;

  /**
   * Carries out the correction (post-force) integration of the motion of all
   * particles by using the acceleration with the explicit Euler method.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param momentum Momentum of particles
   * @param force Force acting on particles
   * @param MOI A container of moment of inertia of particles
   * @param acceleration A container of acceleration of particles
   */
  virtual void
  integrate_post_force(
    Particles::ParticleHandler<dim> &        particle_handler,
    Tensor<1, dim>                           body_force,
    double                                   time_step,
    std::unordered_map<int, Tensor<1, dim>> &momentum,
    std::unordered_map<int, Tensor<1, dim>> &force,
    std::unordered_map<int, double> &        MOI,
    std::unordered_map<int, Tensor<1, dim>> &acceleration) override;
};

#endif
