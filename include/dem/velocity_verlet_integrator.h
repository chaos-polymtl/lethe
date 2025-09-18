// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_velocity_verlet_integrator_h
#define lethe_velocity_verlet_integrator_h

#include <dem/dem_solver_parameters.h>
#include <dem/integrator.h>

using namespace dealii;

/**
 * @brief Implementation of a classical velocity verlet scheme for the
 * integration of the particle motion. Note that reinitialization of force and
 * torque is also integrated into integration class.
 *
 * Velocity Verlet is a second-order integration scheme.
 *
 * Calculation procedure:
 *
 * Calculation of half-step velocity
 * v(n+1/2) = v(n)     + 1/2 * a(n)   * dt
 *
 * Updating particle position
 * x(n+1)   = x(n)     + v(n+1/2)     * dt
 *
 * Updating acceleration and velocity of particle
 * a(n+1)   = F(n+1) / m
 * v(n+1)   = v(n+1/2) + 1/2 * a(n+1) * dt
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 */
template <int dim, typename PropertiesIndex>
class VelocityVerletIntegrator : public Integrator<dim, PropertiesIndex>
{
public:
  VelocityVerletIntegrator()
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
   * @brief Calculate the integration of the motion of all
   * particles by using the Velocity Verlet method.
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

  void
  integrate_with_advected_particles(
    Particles::ParticleHandler<dim>                 &particle_handler,
    const Tensor<1, 3>                              &body_force,
    const double                                     time_step,
    std::vector<Tensor<1, 3>>                       &torque,
    std::vector<Tensor<1, 3>>                       &force,
    const std::vector<double>                       &MOI,
    const parallel::distributed::Triangulation<dim> &triangulation,
    AdaptiveSparseContacts<dim, PropertiesIndex>    &sparse_contacts_object);
};

#endif
