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

#include <dem/dem_solver_parameters.h>
#include <dem/pw_contact_info_struct.h>
#include <dem/pw_contact_force.h>

#include <deal.II/distributed/tria.h>

using namespace dealii;

#ifndef grid_motion_h
#  define grid_motion_h

/**
 * Carries out motion of the triangulation, including rotational, translational
 * and rotational-translational motions.
 *
 * @note Grid motion with a ParsedFunction should be added to the code as a future
 * improvement
 * @param triangulation Triagulation
 * @param dem_parameters Input DEM parameters in the parameter handler file
 *
 */

template <int dim>
class GridMotion
{
  using FuncPtrType = void (GridMotion<dim>::*)(
    parallel::distributed::Triangulation<dim> &triangulation);
  FuncPtrType grid_motion;

public:
  /**
   * The constructor sets up the grid motion function defined in the parameter
   * handler, and calculates the rotation angle for rotational and
   * translational-rotational motions.
   *
   * @param dem_parameters DEM parameters defined using the parameter handler
   * @param dem_time_step DEM time-step
   */
  GridMotion(const DEMSolverParameters<dim> &dem_parameters,
             const double &                  dem_time_step,
             std::shared_ptr<PWContactForce<dim>> pw_contact_force_object={});

  /**
   * Calls the desired grid motion.
   *
   * @param triangulation Triangulation
   */
  void
  move_grid(parallel::distributed::Triangulation<dim> &triangulation)
  {
    (this->*grid_motion)(triangulation);
  }

  /**
   * Carries out updating the boundary points and normal vectors in the
   * particle-wall contact list.
   *
   * @param pw_pairs_in_contact The particle-wall contact list container.
   * We will update the positions of the boundary points and normal vectors
   * directly in this container.
   * @param updated_boundary_points_and_normal_vectors A map that contains
   * updated points on boundaries and normal vectors of the boundary faces.
   * This container is used when the grid is moving. We use this vector to
   * update the boundary points and normal vectors in particle-wall contact
   * list.
   */
  void
  update_boundary_points_and_normal_vectors_in_contact_list(
    std::unordered_map<
      types::particle_index,
      std::map<types::particle_index, pw_contact_info_struct<dim>>>
      &pw_pairs_in_contact,
    const std::map<unsigned int, std::pair<Tensor<1, dim>, Point<dim>>>
      &updated_boundary_points_and_normal_vectors);

private:
  /**
   * Carries out rotational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_rotational(
    parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * Carries out translational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_translational(
    parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * Carries out translational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_due_particles_forces(
    parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * Update member information like triangulation_forces
   * triangulation_forces variable is used to say at this moment t of the simulation,
   * this is the forces applied to the triangulation. But the information,
   * when accessed, is in a form of a map of tensors with as keys the boundary_id.
   * But we only need the total forces to deduce the motion displacement.
   * New update in the future like triangulation_torques would be needed to simulate
   * rotation motion.
   */
  void
  update_parameters_before_motion();

  /**
   * Update member information like inertia tensor,
   * center of mass coordinate.
   * After rotation, the inertia around static axis would change.
   * After translation, the triangulation's center of mass would change,
   * so we need to update them.
   */
  void
  update_parameters_after_motion();

  /**
   * Calculate motion_parameters as rotation_angle and shift_vector
   * members
   */
  void
  calculate_motion_parameters();

  // Since the DEM time-step and rotational speed are constant, we calculate the
  // rotation angle at each time-step once in the constructor and define it as a
  // member variable.
  Tensor<1,3> rotation_angle;
  Tensor<1,dim> shift_vector;

  // Grid motion class need some members to permit the motion
  const double                         triangulation_mass;
  Tensor<1,3>                         triangulation_inertia;
  Tensor<1,dim>                       triangulation_forces;
  Tensor<1,dim>                       triangulation_torques;
  const double dt;
  Tensor<1,dim>     boundary_translational_velocity;
  Tensor<1,3> boundary_rotational_velocity;
  std::shared_ptr<PWContactForce<dim>> pw_contact_force_object;
};


#endif /* grid_motion_h */
