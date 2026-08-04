// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_integrator_h
#define lethe_integrator_h

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

using namespace dealii;

/**
 * @brief Base interface for classes that carry out the integration of the velocity and
 * position of particles with inertia
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class Integrator
{
public:
  /**
   * @brief The constructor to the integrator class is currently blank
   * Eventually it might be a good idea to have the integration class contain
   * the index to the velocity property, the force property and the acceleration
   * property manually
   */
  Integrator()
  {}

  virtual ~Integrator()
  {}

  /**
   * @brief Indicate if the integration scheme stores the velocities staggered by
   * half a time step with respect to the positions. Such schemes are not
   * self-starting: they require an opening half velocity step
   * (integrate_half_step_location) and a closing half velocity step
   * (integrate_half_step_velocity) to synchronize the velocities with the
   * positions at the final time.
   *
   * @return true if the scheme requires an opening and a closing half step.
   */
  virtual bool
  is_half_step_required() const = 0;

  /**
   * @brief Integrate the opening half step of a staggered scheme. The velocities
   * are advanced by half a time step and the positions by a full time step,
   * which brings the velocities to t(1/2) and the positions to t(1).
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
    const std::vector<double>       &MOI) = 0;

  /**
   * @brief Integrate the closing half step of a staggered scheme. The velocities
   * are advanced by half a time step, which synchronizes them with the
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
    const std::vector<double>       &MOI) = 0;

  /**
   * @brief Integrate the closing half step of a staggered scheme when the
   * adaptive sparse contacts are enabled. Only the particles located in mobile
   * cells have their velocities advanced.
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
    AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object) = 0;

  /**
   * @brief Integrate particles' velocity and position.
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
            const std::vector<double>       &MOI) = 0;

  virtual void
  integrate(
    Particles::ParticleHandler<dim>                 &particle_handler,
    const Tensor<1, 3>                              &body_force,
    const double                                     time_step,
    std::vector<Tensor<1, 3>>                       &torque,
    std::vector<Tensor<1, 3>>                       &force,
    const std::vector<double>                       &MOI,
    const parallel::distributed::Triangulation<dim> &triangulation,
    AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object) = 0;
};

#endif
