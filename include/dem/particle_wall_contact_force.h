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
#include <core/auxiliary_math_functions.h>
#include <core/dem_properties.h>
#include <core/serial_solid.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_contact_info.h>

#include <boost/math/special_functions.hpp>
#include <boost/range/adaptor/map.hpp>

#include <cmath>
#include <iostream>

using namespace dealii;

#ifndef particle_wall_contact_force_h
#  define particle_wall_contact_force_h

/**
 * Base interface for classes that carry out the calculation of particle-wall
 * contact force
 */

template <int dim>
class ParticleWallContactForce
{
public:
  ParticleWallContactForce(const DEMSolverParameters<dim> &dem_parameters)
    : n_particle_types(
        dem_parameters.lagrangian_physical_properties.particle_type_number)
  {
    effective_youngs_modulus.resize(n_particle_types);
    effective_shear_modulus.resize(n_particle_types);
    effective_coefficient_of_restitution.resize(n_particle_types);
    effective_coefficient_of_friction.resize(n_particle_types);
    effective_coefficient_of_rolling_friction.resize(n_particle_types);
    effective_surface_energy.resize(n_particle_types);
    model_parameter_beta.resize(n_particle_types);

    this->boundary_translational_velocity_map =
      dem_parameters.boundary_conditions.boundary_translational_velocity;
    this->boundary_rotational_speed_map =
      dem_parameters.boundary_conditions.boundary_rotational_speed;
    this->boundary_rotational_vector =
      dem_parameters.boundary_conditions.boundary_rotational_vector;
    this->point_on_rotation_vector =
      dem_parameters.boundary_conditions.point_on_rotation_axis;
  }

  virtual ~ParticleWallContactForce()
  {}

  /**
   * Carries out the calculation of the particle-wall contact force using the
   * contact pair information obtained from the particle-wall fine search and
   * physical properties of particles and walls
   *
   * @param particle_wall_pairs_in_contact Required information for the calculation of the
   * particle-wall contact force
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   */
  virtual void
  calculate_particle_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_wall_in_contact
                              &particle_wall_pairs_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force) = 0;

  /**
   * Carries out the calculation of particle-floating mesh contact force using
   * the contact pair container
   *
   * @param particle_floating_mesh_in_contact A container that stores the information of
   * particle-floating mesh contact
   * @param dt DEM time step
   * @param torque Torque acting on particles
   * @param force Force acting on particles
   * @param solids Floating solids
   */
  virtual void
  calculate_particle_floating_wall_contact_force(
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_in_contact
                              &particle_floating_mesh_in_contact,
    const double               dt,
    std::vector<Tensor<1, 3>> &torque,
    std::vector<Tensor<1, 3>> &force,
    const std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> &solids) = 0;

  std::map<types::boundary_id, Tensor<1, 3>>
  get_force()
  {
    return force_on_walls;
  }

  std::map<types::boundary_id, Tensor<1, 3>>
  get_torque()
  {
    return torque_on_walls;
  }

  /** This function is used to find the projection of vector_a on
   * vector_b
   * @param vector_a A vector which is going to be projected on vector_b
   * @param vector_b The projection vector of vector_a
   * @return The projection of vector_a on vector_b
   */
  inline Tensor<1, 3>
  find_projection(const Tensor<1, 3> &vector_a, const Tensor<1, 3> &vector_b)
  {
    Tensor<1, 3> vector_c;
    vector_c = ((vector_a * vector_b) / (vector_b.norm_square())) * vector_b;

    return vector_c;
  }

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
    particle_wall_contact_info<dim> &contact_pair_information,
    const Point<3>                  &particle_position,
    const ArrayView<const double>   &particle_properties,
    const double                     dt);

  /**
   * Carries out updating the contact pair information for particle-floating
   * wall contacts
   *
   * @param contact_pair_information Contact information of a particle-wall pair
   * in neighborhood
   * @param particle_properties Properties of particle in contact
   * @param dt DEM time step
   * @param cut_cell_translational_velocity Translational velocity of the cut cell
   * @param cut_cell_rotational_veclocity Rotational veclocity of the cut cell
   * @param center_of_rotation_particle_distance distance between particle and
   * the center of rotation of the floating mesh
   */
  void
  update_particle_floating_wall_contact_information(
    particle_wall_contact_info<dim> &contact_pair_information,
    const ArrayView<const double>   &particle_properties,
    const double                     dt,
    const Tensor<1, 3>              &cut_cell_translational_velocity,
    const Tensor<1, 3>              &cut_cell_rotational_velocity,
    const double                     center_of_rotation_particle_distance);

  /**
   * Carries out applying the calculated force and torque on the particle in
   * contact with the given wall, for both non-linear and linear contact force
   * calculations
   *
   * @param forces_and_torques A tuple which contains: 1, normal force, 2,
   * tangential force, 3, tangential torque and 4, rolling resistance torque of
   * a contact pair
   * @param particle_torque Torque acting on particle
   * @param particle_force Force acting on particle
   * @param point_on_boundary Contact point on the wall
   * @param boundary_id ID of the boundary
   */
  inline void
  apply_force_and_torque(
    const std::tuple<Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>, Tensor<1, 3>>
                   &forces_and_torques,
    Tensor<1, 3>   &particle_torque,
    Tensor<1, 3>   &particle_force,
    const Point<3> &point_on_boundary,
    int             boundary_id = 0)
  {
    // Getting the values from the forces_and_torques tuple, which are: 1,
    // normal force, 2, tangential force, 3, tangential torque and 4, rolling
    // resistance torque
    Tensor<1, 3> normal_force              = std::get<0>(forces_and_torques);
    Tensor<1, 3> tangential_force          = std::get<1>(forces_and_torques);
    Tensor<1, 3> tangential_torque         = std::get<2>(forces_and_torques);
    Tensor<1, 3> rolling_resistance_torque = std::get<3>(forces_and_torques);

    // Calculation of total force
    Tensor<1, 3> total_force = normal_force + tangential_force;

    calculate_force_and_torque_on_boundary(boundary_id,
                                           total_force,
                                           point_on_boundary);

    // Updating the force of particles in the particle handler
    // Since the forces were calculated on the wall, we use the -= operator
    particle_force -= total_force;

    // Updating the torque acting on particles
    // The torque was direcly calculated on the particle, thus we use the +=
    // operator
    particle_torque += tangential_torque + rolling_resistance_torque;
  }

  /** This function is used to calculate the total force and total torque on
   * each boundary
   */
  void
  calculate_force_and_torque_on_boundary(const unsigned int boundary_id,
                                         Tensor<1, 3>       add_force,
                                         const Point<3>     point_contact);

  /** This function is used to initialize a map of vectors to zero
   * with the member class boundary index which has the keys as information
   */
  std::map<unsigned int, Tensor<1, 3>>
  initialize();

  /** This function sums all the forces and torques from all the
   * MPI processes
   */
  void
  mpi_correction_over_calculation_of_forces_and_torques();

  double effective_radius;
  double effective_mass;
  std::unordered_map<unsigned int, Tensor<1, 3>>
                                           boundary_translational_velocity_map;
  std::unordered_map<unsigned int, double> boundary_rotational_speed_map;
  std::unordered_map<unsigned int, Tensor<1, 3>> boundary_rotational_vector;
  std::unordered_map<unsigned int, Point<3>>     point_on_rotation_vector;

  unsigned int n_particle_types;

  std::vector<double> effective_youngs_modulus;
  std::vector<double> effective_shear_modulus;
  std::vector<double> effective_coefficient_of_restitution;
  std::vector<double> effective_coefficient_of_friction;
  std::vector<double> effective_coefficient_of_rolling_friction;
  std::vector<double> effective_surface_energy;
  std::vector<double> model_parameter_beta;

  std::map<unsigned int, Tensor<1, 3>> force_on_walls;
  std::map<unsigned int, Tensor<1, 3>> torque_on_walls;

  bool                            calculate_force_torque_on_boundary;
  Point<3>                        center_mass_container;
  std::vector<types::boundary_id> boundary_index;
  const unsigned int              vertices_per_triangle = 3;
};

#endif /* particle_wall_contact_force_h */
