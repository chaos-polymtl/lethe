// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_explicit_euler_integrator_h
#define lethe_explicit_euler_integrator_h

#include <dem/integrator.h>

using namespace dealii;

/**
 * @brief Implementation of a classical explicit euler scheme for the integration
 * of the particle motion. Note that reinitialization of force and torque is
 * also integrated into integration class
 *
 * @note Euler is a first-order integration scheme. Calculation proceudre:
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * x(n+1) = x(n) + v(n) * dt
 * v(n+1) = v(n) + a(n) * dt
 * a(n+1) = F(n+1) / m
 */
template <int dim, typename PropertiesIndex>
class ExplicitEulerIntegrator : public Integrator<dim, PropertiesIndex>
{
public:
  ExplicitEulerIntegrator()
  {}

  /**
   * @brief Integrate new particles' location after insertion.
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
    const Tensor<1, 3>              &body_force,
    const double                     time_step,
    const std::vector<Tensor<1, 3>> &torque,
    const std::vector<Tensor<1, 3>> &force,
    const std::vector<double>       &MOI) override;

  /**
   * @brief Integrate motion of all particles by using the acceleration with
   * the explicit Euler method.
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
            const Tensor<1, 3>              &body_force,
            const double                     time_step,
            std::vector<Tensor<1, 3>>       &torque,
            std::vector<Tensor<1, 3>>       &force,
            const std::vector<double>       &MOI) override;

  virtual void
  integrate(Particles::ParticleHandler<dim>                 &particle_handler,
            const Tensor<1, 3>                              &body_force,
            const double                                     time_step,
            std::vector<Tensor<1, 3>>                       &torque,
            std::vector<Tensor<1, 3>>                       &force,
            const std::vector<double>                       &MOI,
            const parallel::distributed::Triangulation<dim> &triangulation,
            AdaptiveSparseContacts<dim, PropertiesIndex>
              &sparse_contacts_object) override;
};

#endif
