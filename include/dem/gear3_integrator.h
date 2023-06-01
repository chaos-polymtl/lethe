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

#ifndef gear3_integrator_h
#  define gear3_integrator_h

/**
 * Implementation of a gear3 scheme for the integration
 * of the particle motion. Note that reinitilization of force and torque is also
 * integrated into integration class
 *
 * @note Gear3 is a third-order predictor-correcton integration scheme. In gear3,
 * acceleration, velocity and position of particles are updated using the
 * following method. b is the first derivative of acceleration
 *
 * Prediction:
 * b(n+1,prediction) = b(n)
 * a(n+1,prediction) = a(n) + b(n) * dt
 * v(n+1,prediction) = v(n) + a(n) * dt + 1/2 * b(n) * dt ^ 2
 * x(n+1,prediction) = x(n) + v(n) * dt + 1/2 * a(n) * dt ^ 2 + 1/6 * b(n) * dt
 * ^ 3
 *
 * a(n+1,correction) = F(n+1) / m
 * delata_a(n+1) = a(n+1,correction) - a(n+1,prediction)
 *
 * Correction:
 * ┌                   ┐   ┌                   ┐   ┌                ┐
 * | x(n+1,correction) |   | x(n+1,prediction) |   | eta1 * dt ^ 2  |
 * | v(n+1,correction) | = | v(n+1,prediction) | + | eta2 * dt      | * delata_a(n+1)
 * | a(n+1,correction) |   | a(n+1,prediction) |   | eta3           |
 * | b(n+1,correction) |   | b(n+1,prediction) |   | eta4 * dt ^ -1 |
 * └                   ┘   └                   ┘   └                ┘
 * where eta1 - eta4 are gear3 method coefficients, eta1 = 1/12, eta2 = 5/12,
 * eta3 = 1, eta4 = 1
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class Gear3Integrator : public Integrator<dim>
{
public:
  Gear3Integrator()
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
   * Carries out the integration of the motion of all
   * particles by using the Gear3 method.
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
  Point<dim>     predicted_location;
  Tensor<1, dim> corrected_accereration;
  Tensor<1, dim> acceleration_deviation;
  Point<dim>     corrected_location;
  double         correction;
};

#endif
