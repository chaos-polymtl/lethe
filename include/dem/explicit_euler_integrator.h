// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
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
   * @brief Indicate that the explicit Euler scheme keeps the velocities and the
   * positions synchronized and thus does not require any half step.
   *
   * @return Always false.
   */
  virtual bool
  is_half_step_required() const override
  {
    return false;
  }

  /**
   * @brief Opening half step of a staggered scheme. The explicit Euler scheme is
   * self-starting, so this function must never be called.
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
    std::vector<Tensor<1, 3>>       &torque,
    std::vector<Tensor<1, 3>>       &force,
    const std::vector<double>       &MOI) override;

  /**
   * @brief Closing half step of a staggered scheme. The velocities of the explicit
   * Euler scheme are always synchronized with the positions, so this function
   * must never be called.
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
  integrate_half_step_velocity(
    Particles::ParticleHandler<dim> &particle_handler,
    const Tensor<1, 3>              &body_force,
    const double                     time_step,
    std::vector<Tensor<1, 3>>       &torque,
    std::vector<Tensor<1, 3>>       &force,
    const std::vector<double>       &MOI) override;

  /**
   * @brief Closing half step of a staggered scheme when the adaptive sparse
   * contacts are enabled. The velocities of the explicit Euler scheme are
   * always synchronized with the positions, so this function must never be
   * called.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param MOI A container of moment of inertia of particles
   * @param triangulation The triangulation of the background mesh
   * @param sparse_contacts_object The object that stores the mobility status of
   * the cells
   */
  virtual void
  integrate_half_step_velocity(
    Particles::ParticleHandler<dim>                 &particle_handler,
    const Tensor<1, 3>                              &body_force,
    const double                                     time_step,
    std::vector<Tensor<1, 3>>                       &torque,
    std::vector<Tensor<1, 3>>                       &force,
    const std::vector<double>                       &MOI,
    const parallel::distributed::Triangulation<dim> &triangulation,
    AdaptiveSparseContacts<dim, PropertiesIndex>    &sparse_contacts_object)
    override;

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
