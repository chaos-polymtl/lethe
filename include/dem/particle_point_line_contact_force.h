/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_point_line_contact_info_struct.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_point_line_force_h
#  define particle_point_line_force_h

/**
 * Calculation of the non-linear particle-point and particle-line contact force
 * using the information obtained from the fine search and physical properties
 * of particles and walls. Since the inputs and the calculation method of
 * particle-point and particle-line contact forces are similar, we only used one
 * function for both tasks
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticlePointLineForce
{
public:
  ParticlePointLineForce<dim>();

  /**
   * Carries out the calculation of the particle-point contact force using
   * non-linear (Hertzian) model
   *
   * @param particle_point_pairs_in_contact Required information for
   * calculation of the particle-point contact force
   * @param physical_properties DEM physical_properties declared in the .prm
   * file
   * @param force Force acting on particles
   */
  void
  calculate_particle_point_contact_force(
    const typename dem_data_containers::dem_data_structures<dim>::
      particle_point_line_contact_info *particle_point_line_pairs_in_contact,
    const Parameters::Lagrangian::LagrangianPhysicalProperties
      &                        lagrangian_physical_properties,
    std::vector<Tensor<1, 3>> &force);

  /**
   * Carries out the calculation of the particle-line contact force using
   * non-linear (Hertzian) model
   *
   * @param particle_line_pairs_in_contact Required information for
   * calculation of the particle-line contact force
   * @param physical_properties DEM physical_properties declared in the .prm
   * file
   * @param force Force acting on particles
   */
  void
  calculate_particle_line_contact_force(
    const typename dem_data_containers::dem_data_structures<
      dim>::particle_point_line_contact_info *particle_line_pairs_in_contact,
    const Parameters::Lagrangian::LagrangianPhysicalProperties
      &                        lagrangian_physical_properties,
    std::vector<Tensor<1, 3>> &force);

private:
  /** This private function is used to find the projection of point_p on
   * a line with beginning and ending vertices of point_a and point_b,
   * respectively
   * @param point_p A point which is going to be projected on the line
   * @param point_a Beginning point of the line
   * @param point_b Ending point of the line
   * @return The projection of point_p on the line (from point_a to point_b)
   */

  Point<3>
  find_projection_point(const Point<3> &point_p,
                        const Point<3> &point_a,
                        const Point<3> &point_b);
};

#endif
