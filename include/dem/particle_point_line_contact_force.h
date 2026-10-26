// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_point_line_force_h
#define lethe_particle_point_line_force_h

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

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

  /**
   * @brief Calculate the displacements that separate the particles from the
   * boundary points they overlap. This "contact force" model is used by the
   * packed insertion method.
   *
   * @param[in] particle_point_pairs_in_contact Required information for
   * calculation of the particle-point contacts.
   * @param[in,out] displacement Displacement of each local particle. The
   * displacement of each particle-point contact is added to the entry of its
   * particle.
   */
  void
  calculate_particle_point_shift_displacement(
    const typename DEM::dem_data_structures<dim>::particle_point_in_contact
                              &particle_point_pairs_in_contact,
    std::vector<Tensor<1, 3>> &displacement);

  /**
   * @brief Calculate the displacements that separate the particles from the
   * boundary lines they overlap. This "contact force" model is used by the
   * packed insertion method.
   *
   * @param[in] particle_line_pairs_in_contact Required information for
   * calculation of the particle-line contacts.
   * @param[in,out] displacement Displacement of each local particle. The
   * displacement of each particle-line contact is added to the entry of its
   * particle.
   */
  void
  calculate_particle_line_shift_displacement(
    const typename DEM::dem_data_structures<dim>::particle_line_in_contact
                              &particle_line_pairs_in_contact,
    std::vector<Tensor<1, 3>> &displacement);

  /**
   * @brief Return the number of particle-point and particle-line contacts
   * that occurred in the present pseudo-time step when using the packed
   * insertion method.
   */
  unsigned int
  get_number_of_contacts() const
  {
    return n_contacts;
  }

  /**
   * @brief Reset the counter of number of contacts for the next pseudo-time step.
   */
  void
  reset_number_of_contacts()
  {
    n_contacts = 0;
  }

private:
  /**
   * @brief Add the displacement that separates a particle from the boundary
   * point or line it overlaps to the displacement of the particle. The
   * displacement is along the normal from the closest point of the boundary
   * to the particle center and is slightly larger than the normal overlap, as
   * in the particle-wall shift model.
   *
   * @param[in] particle Iterator of the particle in contact.
   * @param[in] particle_location Location of the particle center in 3D.
   * @param[in] contact_point Point of the boundary point or line that is the
   * closest to the particle center.
   * @param[in,out] displacement Displacement of each local particle.
   */
  void
  add_shift_displacement(const Particles::ParticleIterator<dim> &particle,
                         const Point<3>            &particle_location,
                         const Point<3>            &contact_point,
                         std::vector<Tensor<1, 3>> &displacement);

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

  /**
   * @brief Number of particle-point and particle-line contacts in the present
   * pseudo-time step when using the packed insertion method.
   */
  unsigned int n_contacts = 0;
};

#endif
