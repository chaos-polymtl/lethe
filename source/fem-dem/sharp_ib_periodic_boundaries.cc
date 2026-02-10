// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <fem-dem/sharp_ib_periodic_boundaries.h>

#include <limits>

template <int dim>
SharpIBPeriodicBoundaries<dim>::SharpIBPeriodicBoundaries()
  : periodic_enabled(false)
  , dem_periodic_enabled(false)
  , cfd_periodic_enabled(false)
  , periodicity_mismatch_detected(false)
  , periodic_boundary_0(0)
  , periodic_boundary_1(0)
  , periodic_direction(0)
  , periodic_offset_magnitude(0.0)
{}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::setup(
  const bool                                    dem_periodic_enabled_input,
  const types::boundary_id                      dem_periodic_boundary_0_input,
  const types::boundary_id                      dem_periodic_boundary_1_input,
  const unsigned int                            dem_periodic_direction_input,
  const BoundaryConditions::BoundaryConditions &cfd_boundary_conditions,
  const ConditionalOStream                     &pcout)
{
  dem_periodic_enabled = dem_periodic_enabled_input;

  // Check if CFD has periodic boundaries
  cfd_periodic_enabled = false;
  types::boundary_id cfd_periodic_boundary_0 = 0;
  types::boundary_id cfd_periodic_boundary_1 = 0;
  unsigned int       cfd_periodic_direction  = 0;

  for (auto const &[id, type] : cfd_boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          cfd_periodic_enabled      = true;
          cfd_periodic_boundary_0   = id;
          cfd_periodic_boundary_1   = cfd_boundary_conditions.periodic_neighbor_id.at(id);
          cfd_periodic_direction    = cfd_boundary_conditions.periodic_direction.at(id);
          break; // Currently support only one periodic direction
        }
    }

  // Determine if we can enable periodicity for sharp IB coupling
  periodic_enabled = dem_periodic_enabled || cfd_periodic_enabled;

  if (!periodic_enabled)
    {
      // No periodicity configured anywhere
      return;
    }

  // Check consistency
  if (dem_periodic_enabled && cfd_periodic_enabled)
    {
      // Both are enabled - check if they match
      if (dem_periodic_boundary_0_input != cfd_periodic_boundary_0 ||
          dem_periodic_boundary_1_input != cfd_periodic_boundary_1 ||
          dem_periodic_direction_input != cfd_periodic_direction)
        {
          periodicity_mismatch_detected = true;
        }
      else
        {
          periodicity_mismatch_detected = false;
        }

      // Use CFD configuration as the source of truth for coupling
      // (since refinement is on the fluid mesh)
      periodic_boundary_0 = cfd_periodic_boundary_0;
      periodic_boundary_1 = cfd_periodic_boundary_1;
      periodic_direction  = cfd_periodic_direction;
    }
  else if (cfd_periodic_enabled)
    {
      // Only CFD is periodic - use its configuration
      periodic_boundary_0           = cfd_periodic_boundary_0;
      periodic_boundary_1           = cfd_periodic_boundary_1;
      periodic_direction            = cfd_periodic_direction;
      periodicity_mismatch_detected = true; // DEM not periodic
    }
  else
    {
      // Only DEM is periodic - use its configuration
      periodic_boundary_0           = dem_periodic_boundary_0_input;
      periodic_boundary_1           = dem_periodic_boundary_1_input;
      periodic_direction            = dem_periodic_direction_input;
      periodicity_mismatch_detected = true; // CFD not periodic
    }

  check_consistency_and_warn(pcout);
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::set_periodic_offset(
  const Tensor<1, dim> &offset)
{
  periodic_offset           = offset;
  periodic_offset_magnitude = offset.norm();
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::check_consistency_and_warn(
  const ConditionalOStream &pcout)
{
  if (!periodic_enabled)
    return;

  if (periodicity_mismatch_detected)
    {
      pcout << "***********************************************" << std::endl;
      pcout << "WARNING: Sharp IB Periodic Boundary Mismatch" << std::endl;
      pcout << "***********************************************" << std::endl;
      pcout << "DEM periodic enabled: " << dem_periodic_enabled << std::endl;
      pcout << "CFD periodic enabled: " << cfd_periodic_enabled << std::endl;

      if (dem_periodic_enabled && cfd_periodic_enabled)
        {
          pcout << "Periodic boundaries are enabled for both DEM and CFD, but "
                   "with different configurations."
                << std::endl;
        }
      else if (dem_periodic_enabled)
        {
          pcout << "Periodic boundaries enabled for DEM but NOT for CFD."
                << std::endl;
        }
      else
        {
          pcout << "Periodic boundaries enabled for CFD but NOT for DEM."
                << std::endl;
        }

      pcout
        << "This mismatch may lead to physically inconsistent coupling behavior."
        << std::endl;
      pcout << "For proper sharp IB coupling with periodicity, ensure both DEM "
               "and CFD"
            << std::endl;
      pcout << "periodic boundaries are configured identically." << std::endl;
      pcout << "***********************************************" << std::endl;
    }
  else
    {
      pcout << "Sharp IB periodic boundaries configured successfully."
            << std::endl;
      pcout << "  Periodic direction: " << periodic_direction << std::endl;
      pcout << "  Periodic boundary IDs: " << periodic_boundary_0 << ", "
            << periodic_boundary_1 << std::endl;
    }
}

template <int dim>
Point<dim>
SharpIBPeriodicBoundaries<dim>::wrap_point(const Point<dim> &point,
                                            const double      domain_lower,
                                            const double      domain_upper) const
{
  if (!periodic_enabled)
    return point;

  Point<dim> wrapped_point = point;

  // Check if point is outside periodic domain
  const double coord = point[periodic_direction];

  if (coord < domain_lower)
    {
      // Wrap from lower to upper
      wrapped_point += periodic_offset;
    }
  else if (coord > domain_upper)
    {
      // Wrap from upper to lower
      wrapped_point -= periodic_offset;
    }

  return wrapped_point;
}

template <int dim>
std::vector<Point<dim>>
SharpIBPeriodicBoundaries<dim>::generate_periodic_image_points(
  const Point<dim> &point,
  const double      support_radius,
  const double      domain_lower,
  const double      domain_upper) const
{
  std::vector<Point<dim>> images;

  if (!periodic_enabled)
    return images;

  const double coord = point[periodic_direction];

  // Check if stencil support crosses lower periodic boundary
  if (coord - support_radius < domain_lower)
    {
      // Need an image on the upper side
      Point<dim> image_upper = point + periodic_offset;
      images.push_back(image_upper);
    }

  // Check if stencil support crosses upper periodic boundary
  if (coord + support_radius > domain_upper)
    {
      // Need an image on the lower side
      Point<dim> image_lower = point - periodic_offset;
      images.push_back(image_lower);
    }

  return images;
}

template <int dim>
bool
SharpIBPeriodicBoundaries<dim>::is_point_near_periodic_boundary(
  const Point<dim> &point,
  const double      tolerance,
  const double      domain_lower,
  const double      domain_upper) const
{
  if (!periodic_enabled)
    return false;

  const double coord = point[periodic_direction];

  return (std::abs(coord - domain_lower) < tolerance ||
          std::abs(coord - domain_upper) < tolerance);
}

template <int dim>
Point<dim>
SharpIBPeriodicBoundaries<dim>::map_dof_across_periodic_boundary(
  const Point<dim> &dof_location) const
{
  // For now, just add the periodic offset
  // The direction is determined by which boundary we're crossing
  // This is a simplified version - more sophisticated mapping may be needed
  return dof_location + periodic_offset;
}

// Explicit template instantiations
template class SharpIBPeriodicBoundaries<2>;
template class SharpIBPeriodicBoundaries<3>;
