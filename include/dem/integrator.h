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
   * @brief Integrate the opening step of the scheme, which is carried out once,
   * at the very first time step. Whatever the scheme, this advances the
   * positions by a full time step, exactly like integrate() does. Schemes that
   * are not self-starting, that is, schemes that store the velocities staggered
   * by half a time step with respect to the positions, use it to bring their
   * velocities to t(1/2) while the positions reach t(1). Self-starting schemes
   * simply carry out a regular step.
   *
   * Calling this at the first time step is unconditional: it is up to the
   * scheme, and not to the solver, to decide what an opening step means.
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
  integrate_start(Particles::ParticleHandler<dim> &particle_handler,
                  const Tensor<1, 3>              &body_force,
                  const double                     time_step,
                  std::vector<Tensor<1, 3>>       &torque,
                  std::vector<Tensor<1, 3>>       &force,
                  const std::vector<double>       &MOI) = 0;

  /**
   * @brief Integrate the closing step of the scheme, which is carried out once the
   * final position is reached. The positions are never modified. Schemes that
   * store the velocities staggered by half a time step advance them by half a
   * time step here, which synchronizes them with the positions. Schemes whose
   * velocities are already synchronized leave them untouched.
   *
   * Whatever the scheme, the force and the torque containers are reset, since
   * they have been filled by the force evaluation that precedes this call.
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
  integrate_end(Particles::ParticleHandler<dim> &particle_handler,
                const Tensor<1, 3>              &body_force,
                const double                     time_step,
                std::vector<Tensor<1, 3>>       &torque,
                std::vector<Tensor<1, 3>>       &force,
                const std::vector<double>       &MOI) = 0;

  /**
   * @brief Integrate the closing step of the scheme when the adaptive sparse
   * contacts are enabled. Only the particles located in mobile cells have their
   * velocities advanced, the other ones only have their force and torque reset.
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
  integrate_end(
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
