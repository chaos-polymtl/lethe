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

#include <deal.II/distributed/tria.h>

using namespace dealii;

#ifndef grid_motion_h
#  define grid_motion_h

/**
 * Carries out motion of the triangulation, including rotational, translational
 * and rotational-translational motions.
 *
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
             const double &                  dem_time_step);

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
   * Carries out translational-rotational motion of the triangulation
   *
   * @param triangulation Triangulation
   */
  void
  move_grid_translational_rotational(
    parallel::distributed::Triangulation<dim> &triangulation);

  // Since the DEM time-step and rotational speed are constant, we calculate the
  // rotation angle at each time-step once in the constructor and define it as a
  // member variable.
  double rotation_angle;

  unsigned int rotation_axis;

  Tensor<1, dim> shift_vector;
};


#endif /* grid_motion_h */
