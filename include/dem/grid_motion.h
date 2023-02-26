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

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_contact_info.h>

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

template <int dim, int spacedim = dim>
class GridMotion
{
  using FuncPtrType = void (GridMotion<dim, spacedim>::*)(
    Triangulation<dim, spacedim> &triangulation);
  FuncPtrType grid_motion;

public:
  /**
   * The constructor sets up the grid motion function defined in the parameter
   * handler, and calculates the rotation angle for rotational and
   * translational-rotational motions.
   *
   * @param grid_motion_parameters DEM parameters defined using the parameter handler
   * @param dem_time_step DEM time-step
   */
  GridMotion(
    const Parameters::Lagrangian::GridMotion<spacedim> &grid_motion_parameters,
    const double                                        dem_time_step);

  /**
   * Calls the desired grid motion.
   *
   * @param triangulation Triangulation
   */
  void
  move_grid(Triangulation<dim, spacedim> &triangulation)
  {
    (this->*grid_motion)(triangulation);
  }

  /**
   * Carries out updating the boundary points and normal vectors in the
   * particle-wall contact list.
   *
   * @param particle_wall_pairs_in_contact The particle-wall contact list container.
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
    typename DEM::dem_data_structures<spacedim>::particle_wall_in_contact
      &particle_wall_pairs_in_contact,
    const typename DEM::dem_data_structures<
      spacedim>::boundary_points_and_normal_vectors
      &updated_boundary_points_and_normal_vectors);

private:
  /**
   * Carries out rotational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_rotational(Triangulation<dim, spacedim> &triangulation);

  /**
   * Carries out translational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_translational(Triangulation<dim, spacedim> &triangulation);

  // Since the DEM time-step and rotational speed are constant, we calculate the
  // rotation angle at each time-step once in the constructor and define it as a
  // member variable.
  double rotation_angle;

  unsigned int rotation_axis;

  Tensor<1, spacedim> shift_vector;
};


#endif /* grid_motion_h */
