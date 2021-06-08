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
#include <boost/math/special_functions.hpp>
#include <boost/range/adaptor/map.hpp>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_info_struct.h>
#include <math.h>

#include <iostream>

using namespace dealii;

#ifndef particle_wall_contact_force_h
#  define particle_wall_contact_force_h

/**
 * Base interface for classes that carry out the calculation of particle-wall
 * contact force
 */

template <int dim>
class PWContactForce
{
public:
  PWContactForce()
  {}

  virtual ~PWContactForce()
  {}

  /**
   * Carries out the calculation of the particle-wall contact force using the
   * contact pair information obtained from the particle-wall fine search and
   * physical properties of particles and walls
   *
   * @param pw_pairs_in_contact Required information for the calculation of the
   * particle-wall contact force
   * @param dt DEM time step
   * @param momentum An unordered_map of momentum of particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_pw_contact_force(
    std::unordered_map<
      types::particle_index,
      std::map<types::particle_index, pw_contact_info_struct<dim>>>
      &           pw_pairs_in_contact,
    const double &dt,
    std::unordered_map<types::particle_index, Tensor<1, dim>> &momentum,
    std::unordered_map<types::particle_index, Tensor<1, dim>> &force) = 0;

  void
  calculate_pw_force_torque(
    std::unordered_map<
      types::particle_index,
      std::map<types::particle_index, pw_contact_info_struct<dim>>>
      &                      pw_pairs_in_contact,
    const double &           dt,
    DEMSolverParameters<dim> parameters);

  Tensor<1, dim> calculation_total_torque(Tensor<1, dim> total_force,
                                          Point<dim>     center_mass,
                                          Point<dim>     point_on_boundary);

protected:
  /**
   * Carries out updating the contact pair information for both non-linear and
   * linear contact force calculations
   *
   * @param contact_pair_information Contact information of a particle-wall pair
   * in neighborhood
   * @param particle_properties Properties of particle in contact
   * @param dt DEM time step
   */
  void
  update_contact_information(
    pw_contact_info_struct<dim> &  contact_pair_information,
    const ArrayView<const double> &particle_properties,
    const double &                 dt);

  /**
   * Carries out applying the calculated force and torque on the local-local
   * particle pair in contact, for both non-linear and linear contact force
   * calculations
   *
   * @param particle_properties Properties of particle in contact with wall
   * @param forces_and_torques A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair
   * @param particle_momentum Momentum of particle
   * @param particle_force Force acting on particle
   */
  inline void
  apply_force_and_torque(const std::tuple<Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>> &forces_and_torques,
                         Tensor<1, dim> &                  particle_momentum,
                         Tensor<1, dim> &                  particle_force)
  {
    // Getting the values from the forces_and_torques tuple, which are: 1,
    // normal force, 2, tangential force, 3, tangential torque and 4, rolling
    // resistance torque
    Tensor<1, dim> normal_force              = std::get<0>(forces_and_torques);
    Tensor<1, dim> tangential_force          = std::get<1>(forces_and_torques);
    Tensor<1, dim> tangential_torque         = std::get<2>(forces_and_torques);
    Tensor<1, dim> rolling_resistance_torque = std::get<3>(forces_and_torques);

    // Calculation of total force
    Tensor<1, dim> total_force = normal_force + tangential_force;

    // Updating the force of particles in the particle handler
    for (int d = 0; d < dim; ++d)
      {
        particle_force[d] = particle_force[d] + total_force[d];
      }

    // Updating the torque acting on particles
    for (int d = 0; d < dim; ++d)
      {
        particle_momentum[d] = particle_momentum[d] + tangential_torque[d] +
                               rolling_resistance_torque[d];
      }
  }

  /** This function is used to find the projection of vector_a on
   * vector_b
   * @param vector_a A vector which is going to be projected on vector_b
   * @param vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  inline Tensor<1, dim>
  find_projection(const Tensor<1, dim> &vector_a,
                  const Tensor<1, dim> &vector_b)
  {
    Tensor<1, dim> vector_c;
    vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

    return vector_c;
  }

  double triangulation_radius;
  double effective_radius;
  double effective_mass;
  std::unordered_map<types::particle_index, Tensor<1, dim>>
    boundary_translational_velocity_map;
  std::unordered_map<types::particle_index, double>
    boundary_rotational_speed_map;
  std::unordered_map<types::particle_index, Tensor<1, dim>>
                                          boundary_rotational_vector;
  std::map<types::particle_index, double> effective_youngs_modulus;
  std::map<types::particle_index, double> effective_shear_modulus;
  std::map<types::particle_index, double> effective_coefficient_of_restitution;
  std::map<types::particle_index, double> effective_coefficient_of_friction;
  std::map<types::particle_index, double>
    effective_coefficient_of_rolling_friction;
};

#endif /* particle_wall_contact_force_h */
