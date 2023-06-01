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

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/disable_contacts.h>

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

#ifndef integration_h
#  define integration_h

/**
 * Base interface for classes that carry out the integration of the velocity and
 * position of particles with inertia
 */

template <int dim>
class Integrator
{
public:
  /**
   * The constructor to the integrator class is currently blank
   * Eventually it might be a good idea to have the integration class contain
   * the index to the velocity property, the force property and the acceleration
   * property manually
   */
  Integrator<dim>()
  {}

  virtual ~Integrator()
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
    const std::vector<double> &      MOI) = 0;

  /**
   * Carries out integrating of particles' velocity and position.
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
            const std::vector<double> &      MOI) = 0;

  virtual void
  integrate(Particles::ParticleHandler<dim> &                particle_handler,
            const Tensor<1, 3> &                             body_force,
            const double                                     time_step,
            std::vector<Tensor<1, 3>> &                      torque,
            std::vector<Tensor<1, 3>> &                      force,
            const std::vector<double> &                      MOI,
            const parallel::distributed::Triangulation<dim> &triangulation,
            DisableContacts<dim> &disable_contacts_object) = 0;
};

#endif /* integration_h */
