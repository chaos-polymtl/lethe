// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_point_line_force_h
#define lethe_particle_point_line_force_h

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <vector>

using namespace dealii;

/**
 * @brief Calculate the non-linear particle-point and particle-line contact force
 * using the information obtained from the fine search and physical properties
 * of particles and walls. Since the inputs and the calculation method of
 * particle-point and particle-line contact forces are similar, we only used one
 * function for both tasks.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class ParticlePointLineForce
{
public:
  ParticlePointLineForce();

  /**
   * @brief Calculate of the particle-point contact force using
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
    const typename DEM::dem_data_structures<dim>::particle_point_in_contact
      *particle_point_pairs_in_contact,
    const Parameters::Lagrangian::LagrangianPhysicalProperties
                              &physical_properties,
    std::vector<Tensor<1, 3>> &force);

  /**
   * @brief Calculate of the particle-line contact force using
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
    const typename DEM::dem_data_structures<dim>::particle_line_in_contact
      *particle_line_pairs_in_contact,
    const Parameters::Lagrangian::LagrangianPhysicalProperties
                              &physical_properties,
    std::vector<Tensor<1, 3>> &force);

private:
  /**
   * @brief This private function is used to find the projection of point_p on
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
