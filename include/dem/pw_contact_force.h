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


#ifndef particle_wall_contact_force_h
#define particle_wall_contact_force_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/tensor.h>

#include <dem/pw_contact_info_struct.h>

#include <unordered_map>

using namespace dealii;

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
   */
  virtual void
  calculate_pw_contact_force(
    std::unordered_map<int, std::map<int, pw_contact_info_struct<dim>>>
      &           pw_pairs_in_contact,
    const double &dt) = 0;

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
   */
  void
  apply_force_and_torque(ArrayView<double> &               particle_properties,
                         const std::tuple<Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>,
                                          Tensor<1, dim>> &forces_and_torques);

  /** This function is used to find the projection of vector_a on
   * vector_b
   * @param vector_a A vector which is going to be projected on vector_b
   * @param vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  Tensor<1, dim>
  find_projection(const Tensor<1, dim> &vector_a,
                  const Tensor<1, dim> &vector_b);

  double                                  triangulation_radius;
  double                                  effective_radius;
  double                                  effective_mass;
  std::unordered_map<int, Tensor<1, dim>> boundary_translational_velocity_map;
  std::unordered_map<int, double>         boundary_rotational_speed_map;
  std::unordered_map<int, Tensor<1, dim>> boundary_rotational_vector;
  std::map<int, double>                   effective_youngs_modulus;
  std::map<int, double>                   effective_shear_modulus;
  std::map<int, double>                   effective_coefficient_of_restitution;
  std::map<int, double>                   effective_coefficient_of_friction;
  std::map<int, double> effective_coefficient_of_rolling_friction;
};

#endif /* particle_wall_contact_force_h */
