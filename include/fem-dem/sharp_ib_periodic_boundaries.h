// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_sharp_ib_periodic_boundaries_h
#define lethe_sharp_ib_periodic_boundaries_h

#include <core/boundary_conditions.h>
#include <core/parameters_lagrangian.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/tensor.h>

/**
 * @brief Combined periodic boundary handler for sharp IB coupling.
 *
 * This class validates the periodic configuration required by Sharp-IB and
 * exposes the minimal geometric helpers needed once that configuration has
 * been accepted.
 *
 * @tparam dim Spatial dimension
 */
template <int dim>
class SharpIBPeriodicBoundaries
{
public:
  /**
   * @brief Validated periodic configuration for Sharp-IB coupling.
   *
   * Sharp currently supports at most one periodic boundary pair. The two
   * boundaries are treated as an ordered pair whose direction is the normal
   * direction used for wrapping and image generation.
   */
  struct PeriodicConfiguration
  {
    bool                       enabled    = false;
    dealii::types::boundary_id boundary_0 = 0;
    dealii::types::boundary_id boundary_1 = 0;
    unsigned int               direction  = 0;
  };

  /**
   * @brief Default constructor. Periodicity is disabled until configured.
   */
  SharpIBPeriodicBoundaries();

  /**
   * @brief Validate and store the periodic configuration used by Sharp-IB.
   *
   * Sharp-IB periodic coupling only supports one periodic boundary pair and
   * requires the DEM and CFD configurations to match exactly. Unsupported or
   * inconsistent configurations throw immediately so the solver never enters a
   * partially periodic state.
   *
   * @param dem_boundary_conditions DEM particle boundary conditions.
   * @param cfd_boundary_conditions CFD boundary conditions object
   * @param pcout Parallel output stream used to report the accepted
   * configuration.
   */
  void
  initialize(
    const Parameters::Lagrangian::BCDEM          &dem_boundary_conditions,
    const BoundaryConditions::BoundaryConditions &cfd_boundary_conditions,
    const dealii::ConditionalOStream             &pcout);

  /**
   * @brief Set the periodic offset after triangulation setup.
   *
   * Must be called after the domain bounds of the active triangulation have
   * been computed.
   *
   * @param offset Periodic offset tensor (displacement from boundary 0 to 1)
   */
  void
  set_periodic_offset(const dealii::Tensor<1, dim> &offset);

  /**
   * @brief Check if periodic boundaries are enabled for sharp IB coupling.
   *
   * @return True if periodicity is enabled and consistent enough to use
   */
  inline bool
  is_periodic_enabled() const
  {
    return configuration.enabled;
  }

  /**
   * @brief Get the validated periodic configuration.
   *
   * @return Active periodic configuration.
   */
  inline const PeriodicConfiguration &
  get_configuration() const
  {
    return configuration;
  }

  /**
   * @brief Get the periodic offset tensor.
   *
   * @return Offset vector from periodic boundary 0 to boundary 1.
   */
  inline dealii::Tensor<1, dim>
  get_periodic_offset() const
  {
    return periodic_offset;
  }

  /**
   * @brief Get the periodic boundary IDs.
   *
   * @return Pair of (boundary_0_id, boundary_1_id).
   */
  inline std::pair<dealii::types::boundary_id, dealii::types::boundary_id>
  get_periodic_boundary_ids() const
  {
    return {configuration.boundary_0, configuration.boundary_1};
  }

  /**
   * @brief Get the periodic direction (0=x, 1=y, 2=z).
   *
   * @return Periodic direction index.
   */
  inline unsigned int
  get_periodic_direction() const
  {
    return configuration.direction;
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
  dealii::Point<dim>
  wrap_point(const dealii::Point<dim> &point,
             const double              domain_lower,
             const double              domain_upper) const;

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
  std::vector<dealii::Point<dim>>
  generate_periodic_image_points(const dealii::Point<dim> &point,
                                 const double              support_radius,
                                 const double              domain_lower,
                                 const double              domain_upper) const;

private:
  /**
   * @brief Read the DEM periodic configuration.
   *
   * @param dem_boundary_conditions DEM boundary conditions.
   * @return Parsed DEM periodic configuration.
   */
  PeriodicConfiguration
  parse_dem_periodic_configuration(
    const Parameters::Lagrangian::BCDEM &dem_boundary_conditions) const;

  /**
   * @brief Read the CFD periodic configuration.
   *
   * @param cfd_boundary_conditions CFD boundary conditions.
   * @return Parsed CFD periodic configuration.
   */
  PeriodicConfiguration
  parse_cfd_periodic_configuration(const BoundaryConditions::BoundaryConditions
                                     &cfd_boundary_conditions) const;

  /**
   * @brief Validate that the DEM and CFD periodic configurations are both
   * supported and identical.
   *
   * @param dem_configuration Parsed DEM periodic configuration.
   * @param cfd_configuration Parsed CFD periodic configuration.
   */
  void
  validate_matching_configurations(
    const PeriodicConfiguration &dem_configuration,
    const PeriodicConfiguration &cfd_configuration) const;

  /**
   * @brief Report the validated periodic configuration used by Sharp-IB.
   *
   * @param pcout Parallel output stream for informational output.
   */
  void
  report_configuration(const dealii::ConditionalOStream &pcout) const;

  // Validated periodic configuration.
  PeriodicConfiguration configuration;

  // Periodic offset vector: the displacement from the lower periodic face to
  // the upper periodic face along the periodic direction. Applied as a
  // translation to map points or particles across the periodic boundary.
  dealii::Tensor<1, dim> periodic_offset;
};

#endif
