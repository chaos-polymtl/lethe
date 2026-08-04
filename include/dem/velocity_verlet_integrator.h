// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
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
 * Velocity Verlet is a second-order integration scheme. It is implemented in
 * its staggered (leapfrog) form: the velocities are stored at the half steps
 * and the positions at the integer steps. Consequently, carrying out n time
 * steps requires n + 1 integration operations.
 *
 * Calculation procedure:
 *
 * Opening half step (integrate_half_step_location), carried out once at the
 * first time step
 * v(1/2)     = v(0)       + 1/2 * a(0)   * dt
 * x(1)       = x(0)       + v(1/2)       * dt
 *
 * Regular steps (integrate), carried out at every subsequent time step
 * a(n)       = F(n) / m
 * v(n+1/2)   = v(n-1/2)   + a(n)         * dt
 * x(n+1)     = x(n)       + v(n+1/2)     * dt
 *
 * Closing half step (integrate_half_step_velocity), carried out once the final
 * position is reached to synchronize the velocities with the positions
 * a(n)       = F(n) / m
 * v(n)       = v(n-1/2)   + 1/2 * a(n)   * dt
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
   * @brief Indicate that the Velocity Verlet scheme stores the velocities
   * staggered by half a time step and thus requires an opening and a closing
   * half step.
   *
   * @return Always true.
   */
  virtual bool
  is_half_step_required() const override
  {
    return true;
  }

  /**
   * @brief Integrate the opening half step of the Velocity Verlet scheme. The
   * velocities are advanced by half a time step and the positions by a full
   * time step.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles. It is reset to zero at the end of
   * the integration
   * @param force Force acting on particles. It is reset to zero at the end of
   * the integration
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
   * @brief Integrate the closing half step of the Velocity Verlet scheme. The
   * velocities are advanced by half a time step to synchronize them with the
   * positions. The positions are left untouched.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles. It is reset to zero at the end of
   * the integration
   * @param force Force acting on particles. It is reset to zero at the end of
   * the integration
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
   * @brief Integrate the closing half step of the Velocity Verlet scheme when the
   * adaptive sparse contacts are enabled. Only the particles located in mobile
   * cells have their velocities advanced, the other ones only have their force
   * and torque reset.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   * @param torque Torque acting on particles. It is reset to zero at the end of
   * the integration
   * @param force Force acting on particles. It is reset to zero at the end of
   * the integration
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
