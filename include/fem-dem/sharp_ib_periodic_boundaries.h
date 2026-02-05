// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sharp_ib_periodic_boundaries_h
#define lethe_sharp_ib_periodic_boundaries_h

#include <core/boundary_conditions.h>

#include <dem/periodic_boundaries_manipulator.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

using namespace dealii;

/**
 * @brief Combined periodic boundary handler for sharp IB coupling.
 *
 * This class provides a unified interface for periodic boundary conditions
 * that is aware of both DEM particle periodicity and CFD field periodicity.
 * It ensures consistent periodic wrapping/mapping across the fluid-particle
 * coupling in the sharp immersed boundary solver.
 *
 * Key responsibilities:
 * - Detect and warn about DEM/CFD periodicity mismatches
 * - Provide periodic offset tensor for coupling operations
 * - Generate periodic image points for stencil support
 * - Map points across periodic boundaries consistently
 * - Enable periodic-aware cell/DoF traversal for refinement, interpolation,
 *   and forcing
 *
 * @tparam dim Spatial dimension
 */
template <int dim>
class SharpIBPeriodicBoundaries
{
public:
  /**
   * @brief Default constructor. Periodicity is disabled until configured.
   */
  SharpIBPeriodicBoundaries();

  /**
   * @brief Configure combined periodic boundaries.
   *
   * This function must be called during solver initialization, after both
   * DEM and CFD periodic boundary conditions have been parsed.
   *
   * @param dem_periodic_enabled True if DEM has periodic boundaries enabled
   * @param dem_periodic_boundary_0 DEM periodic boundary ID (first)
   * @param dem_periodic_direction DEM periodic direction (0=x, 1=y, 2=z)
   * @param cfd_boundary_conditions CFD boundary conditions object
   * @param pcout Parallel output stream for warnings
   */
  void
  setup(const bool                                        dem_periodic_enabled,
        const types::boundary_id                          dem_periodic_boundary_0,
        const unsigned int                                dem_periodic_direction,
        const BoundaryConditions::BoundaryConditions     &cfd_boundary_conditions,
        const ConditionalOStream                         &pcout);

  /**
   * @brief Set the periodic offset after triangulation setup.
   *
   * Must be called after periodic cells have been mapped and the offset
   * distance computed (typically from PeriodicBoundariesManipulator).
   *
   * @param offset Periodic offset tensor (displacement from boundary 0 to 1)
   */
  void
  set_periodic_offset(const Tensor<1, dim> &offset);

  /**
   * @brief Check if periodic boundaries are enabled for sharp IB coupling.
   *
   * @return True if periodicity is enabled and consistent enough to use
   */
  inline bool
  is_periodic_enabled() const
  {
    return periodic_enabled;
  }

  /**
   * @brief Get the periodic offset tensor.
   *
   * @return Offset vector from periodic boundary 0 to boundary 1
   */
  inline Tensor<1, dim>
  get_periodic_offset() const
  {
    return periodic_offset;
  }

  /**
   * @brief Get the periodic direction (0=x, 1=y, 2=z).
   *
   * @return Periodic direction index
   */
  inline unsigned int
  get_periodic_direction() const
  {
    return periodic_direction;
  }

  /**
   * @brief Get the periodic boundary IDs.
   *
   * @return Pair of (boundary_0_id, boundary_1_id)
   */
  inline std::pair<types::boundary_id, types::boundary_id>
  get_periodic_boundary_ids() const
  {
    return {periodic_boundary_0, periodic_boundary_1};
  }

  /**
   * @brief Wrap a point to the periodic domain (minimum image).
   *
   * If the point is outside the periodic domain along the periodic direction,
   * wrap it back using the periodic offset.
   *
   * @param point Point to wrap
   * @param domain_lower Lower bound of domain in periodic direction
   * @param domain_upper Upper bound of domain in periodic direction
   * @return Wrapped point
   */
  Point<dim>
  wrap_point(const Point<dim> &point,
             const double      domain_lower,
             const double      domain_upper) const;

  /**
   * @brief Generate periodic image points for a stencil support region.
   *
   * Given a point and a support radius, returns a list of periodic images
   * that may be needed for interpolation/forcing stencils that cross
   * periodic boundaries.
   *
   * @param point Center point of stencil
   * @param support_radius Radius of stencil support
   * @param domain_lower Lower bound of domain in periodic direction
   * @param domain_upper Upper bound of domain in periodic direction
   * @return Vector of image points (may be empty if no images needed)
   */
  std::vector<Point<dim>>
  generate_periodic_image_points(const Point<dim> &point,
                                 const double      support_radius,
                                 const double      domain_lower,
                                 const double      domain_upper) const;

  /**
   * @brief Check if a point is near a periodic boundary (within tolerance).
   *
   * @param point Point to check
   * @param tolerance Distance tolerance for "near"
   * @param domain_lower Lower bound of domain in periodic direction
   * @param domain_upper Upper bound of domain in periodic direction
   * @return True if point is within tolerance of either periodic boundary
   */
  bool
  is_point_near_periodic_boundary(const Point<dim> &point,
                                  const double      tolerance,
                                  const double      domain_lower,
                                  const double      domain_upper) const;

  /**
   * @brief Map a DoF location across the periodic boundary.
   *
   * Given a DoF support point, returns its periodic image if it exists.
   *
   * @param dof_location DoF support point
   * @return Periodic image of the DoF location
   */
  Point<dim>
  map_dof_across_periodic_boundary(const Point<dim> &dof_location) const;

private:
  /**
   * @brief Check consistency between DEM and CFD periodicity and issue
   * warnings.
   *
   * @param pcout Parallel output stream for warnings
   */
  void
  check_consistency_and_warn(const ConditionalOStream &pcout);

  // Configuration flags
  bool periodic_enabled;
  bool dem_periodic_enabled;
  bool cfd_periodic_enabled;
  bool periodicity_mismatch_detected;

  // Periodic boundary information
  types::boundary_id periodic_boundary_0;
  types::boundary_id periodic_boundary_1;
  unsigned int       periodic_direction; // 0=x, 1=y, 2=z

  // Periodic offset (displacement from boundary 0 to boundary 1)
  Tensor<1, dim> periodic_offset;

  // Magnitude of periodic offset
  double periodic_offset_magnitude;
};

#endif
