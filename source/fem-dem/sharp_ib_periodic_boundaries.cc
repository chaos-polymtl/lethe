// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <fem-dem/sharp_ib_periodic_boundaries.h>

#include <algorithm>
#include <sstream>

template <int dim>
SharpIBPeriodicBoundaries<dim>::SharpIBPeriodicBoundaries()
  : configuration()
{}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::initialize(
  const Parameters::Lagrangian::BCDEM          &dem_boundary_conditions,
  const BoundaryConditions::BoundaryConditions &cfd_boundary_conditions,
  const dealii::ConditionalOStream             &pcout)
{
  // Sharp-IB periodic support is intentionally strict: the fluid and particle
  // sides must describe the same single periodic pair. Allowing CFD-only or
  // DEM-only periodicity would leave refinement, wrapping, particle-wall
  // contact, and cut-cell assembly operating in incompatible frames. This was
  // an implicit assumption throughout the original periodic implementation, so
  // we enforce it explicitly here instead of trying to guess which subsystem
  // should define the "real" periodic configuration.
  const PeriodicConfiguration dem_configuration =
    parse_dem_periodic_configuration(dem_boundary_conditions);
  const PeriodicConfiguration cfd_configuration =
    parse_cfd_periodic_configuration(cfd_boundary_conditions);

  validate_matching_configurations(dem_configuration, cfd_configuration);

  if (!dem_configuration.enabled && !cfd_configuration.enabled)
    {
      configuration = {};
      periodic_offset.clear();
      return;
    }

  configuration = dem_configuration;
  report_configuration(pcout);
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::set_periodic_offset(
  const dealii::Tensor<1, dim> &offset)
{
  periodic_offset = offset;
}

template <int dim>
typename SharpIBPeriodicBoundaries<dim>::PeriodicConfiguration
SharpIBPeriodicBoundaries<dim>::parse_dem_periodic_configuration(
  const Parameters::Lagrangian::BCDEM &dem_boundary_conditions) const
{
  PeriodicConfiguration parsed_configuration;

  const unsigned int periodic_boundary_count =
    std::ranges::count(dem_boundary_conditions.bc_types,
                       Parameters::Lagrangian::BCDEM::BoundaryType::periodic);

  AssertThrow(
    periodic_boundary_count <= 1,
    dealii::ExcMessage(
      "Sharp-IB periodic coupling currently supports exactly one DEM periodic "
      "boundary pair."));

  if (periodic_boundary_count == 0)
    return parsed_configuration;

  parsed_configuration.enabled    = true;
  parsed_configuration.boundary_0 = dem_boundary_conditions.periodic_boundary_0;
  parsed_configuration.boundary_1 = dem_boundary_conditions.periodic_boundary_1;
  parsed_configuration.direction  = dem_boundary_conditions.periodic_direction;

  return parsed_configuration;
}

template <int dim>
typename SharpIBPeriodicBoundaries<dim>::PeriodicConfiguration
SharpIBPeriodicBoundaries<dim>::parse_cfd_periodic_configuration(
  const BoundaryConditions::BoundaryConditions &cfd_boundary_conditions) const
{
  PeriodicConfiguration parsed_configuration;

  std::vector<dealii::types::boundary_id> periodic_boundary_ids;
  periodic_boundary_ids.reserve(cfd_boundary_conditions.type.size());

  for (const auto &[boundary_id, boundary_type] : cfd_boundary_conditions.type)
    {
      if (boundary_type == BoundaryConditions::BoundaryType::periodic)
        periodic_boundary_ids.push_back(boundary_id);
    }

  AssertThrow(
    periodic_boundary_ids.size() <= 1,
    dealii::ExcMessage(
      "Sharp-IB periodic coupling currently supports exactly one CFD periodic "
      "boundary pair."));

  if (periodic_boundary_ids.empty())
    return parsed_configuration;

  const auto boundary_id = periodic_boundary_ids.front();

  parsed_configuration.enabled    = true;
  parsed_configuration.boundary_0 = boundary_id;
  parsed_configuration.boundary_1 =
    cfd_boundary_conditions.periodic_neighbor_id.at(boundary_id);
  parsed_configuration.direction =
    cfd_boundary_conditions.periodic_direction.at(boundary_id);

  return parsed_configuration;
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::validate_matching_configurations(
  const PeriodicConfiguration &dem_configuration,
  const PeriodicConfiguration &cfd_configuration) const
{
  if (!dem_configuration.enabled && !cfd_configuration.enabled)
    return;

  std::ostringstream error_message;
  error_message
    << "Sharp-IB periodic coupling requires DEM and CFD periodic boundary "
       "conditions to both be enabled with exactly matching boundary ids and "
       "direction.\n"
    << "DEM periodic: enabled=" << dem_configuration.enabled;

  if (dem_configuration.enabled)
    {
      error_message << ", pair=(" << dem_configuration.boundary_0 << ", "
                    << dem_configuration.boundary_1
                    << "), direction=" << dem_configuration.direction;
    }

  error_message << "\nCFD periodic: enabled=" << cfd_configuration.enabled;
  if (cfd_configuration.enabled)
    {
      error_message << ", pair=(" << cfd_configuration.boundary_0 << ", "
                    << cfd_configuration.boundary_1
                    << "), direction=" << cfd_configuration.direction;
    }

  AssertThrow(dem_configuration.enabled == cfd_configuration.enabled,
              dealii::ExcMessage(error_message.str()));

  AssertThrow(dem_configuration.boundary_0 == cfd_configuration.boundary_0 &&
                dem_configuration.boundary_1 == cfd_configuration.boundary_1 &&
                dem_configuration.direction == cfd_configuration.direction,
              dealii::ExcMessage(error_message.str()));
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::report_configuration(
  const dealii::ConditionalOStream &pcout) const
{
  if (!configuration.enabled)
    return;

  pcout << "Sharp IB periodic boundaries configured successfully." << std::endl;
  pcout << "  Periodic direction: " << configuration.direction << std::endl;
  pcout << "  Periodic boundary IDs: " << configuration.boundary_0 << ", "
        << configuration.boundary_1 << std::endl;
}

template <int dim>
dealii::Point<dim>
SharpIBPeriodicBoundaries<dim>::wrap_point(const dealii::Point<dim> &point,
                                           const double domain_lower,
                                           const double domain_upper) const
{
  if (!configuration.enabled)
    return point;

  dealii::Point<dim> wrapped_point = point;
  const double       coordinate    = point[configuration.direction];

  // Points that strictly overshoot the lower boundary are wrapped upward
  // (toward the upper face), and points that strictly overshoot the upper
  // boundary are wrapped downward. Points that land exactly on either boundary
  // face are left in place because they are still inside the domain.
  if (coordinate < domain_lower)
    wrapped_point += periodic_offset;
  else if (coordinate > domain_upper)
    wrapped_point -= periodic_offset;

  return wrapped_point;
}

template <int dim>
std::vector<dealii::Point<dim>>
SharpIBPeriodicBoundaries<dim>::generate_periodic_image_points(
  const dealii::Point<dim> &point,
  const double              support_radius,
  const double              domain_lower,
  const double              domain_upper) const
{
  std::vector<dealii::Point<dim>> image_points;

  if (!configuration.enabled)
    return image_points;

  const double coordinate = point[configuration.direction];

  if (coordinate - support_radius < domain_lower)
    image_points.push_back(point + periodic_offset);

  if (coordinate + support_radius > domain_upper)
    image_points.push_back(point - periodic_offset);

  return image_points;
}

template class SharpIBPeriodicBoundaries<2>;
template class SharpIBPeriodicBoundaries<3>;
