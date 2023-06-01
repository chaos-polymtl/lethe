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
#include <dem/integrator.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

#ifndef explicit_euler_integrator_h
#  define explicit_euler_integrator_h

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
   * Carries out integrating of new particles' location after insertion.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param MOI A container of moment of inertia of particles
   */
  virtual void
  integrate_half_step_location(
    Particles::ParticleHandler<dim> &particle_handler,
    const Tensor<1, 3> &             body_force,
    const double                     time_step,
    const std::vector<Tensor<1, 3>> &torque,
    const std::vector<Tensor<1, 3>> &force,
    const std::vector<double> &      MOI) override;

  /**
   * Carries out integration of the motion of all
   * particles by using the acceleration with the explicit Euler method.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param MOI A container of moment of inertia of particles
   */
  virtual void
  integrate(Particles::ParticleHandler<dim> &particle_handler,
            const Tensor<1, 3> &             body_force,
            const double                     time_step,
            std::vector<Tensor<1, 3>> &      torque,
            std::vector<Tensor<1, 3>> &      force,
            const std::vector<double> &      MOI) override;

  virtual void
  integrate(Particles::ParticleHandler<dim> &                particle_handler,
            const Tensor<1, 3> &                             body_force,
            const double                                     time_step,
            std::vector<Tensor<1, 3>> &                      torque,
            std::vector<Tensor<1, 3>> &                      force,
            const std::vector<double> &                      MOI,
            const parallel::distributed::Triangulation<dim> &triangulation,
            DisableContacts<dim> &disable_contacts_object) override;

private:
  Tensor<1, 3> acceleration;
};

#endif
