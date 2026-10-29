// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <fem-dem/sharp_ib_periodic_boundaries.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ostream>
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
  // Fluid geometry and particle motion must use the same periodic pair.
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
  AssertThrow(std::isfinite(offset[configuration.direction]) &&
                offset[configuration.direction] > 0.,
              dealii::ExcMessage(
                "The Sharp-IB period must be positive and finite."));
  for (unsigned int d = 0; d < dim; ++d)
    AssertThrow(d == configuration.direction || offset[d] == 0.,
                dealii::ExcMessage(
                  "Sharp-IB requires an axis-aligned translation."));
  periodic_offset = offset;
}

template <int dim>
typename SharpIBPeriodicBoundaries<dim>::PeriodicConfiguration
SharpIBPeriodicBoundaries<dim>::parse_dem_periodic_configuration(
  const Parameters::Lagrangian::BCDEM &dem_boundary_conditions) const
{
  PeriodicConfiguration parsed_configuration;

  AssertThrow(
    dem_boundary_conditions.periodic_boundaries.size() <= 1,
    dealii::ExcMessage(
      "Sharp-IB periodic coupling currently supports exactly one DEM periodic "
      "boundary pair."));

  if (dem_boundary_conditions.periodic_boundaries.empty())
    return parsed_configuration;

  const auto periodic_pair =
    dem_boundary_conditions.periodic_boundaries.begin();
  parsed_configuration.enabled    = true;
  parsed_configuration.boundary_0 = periodic_pair->first;
  parsed_configuration.boundary_1 = periodic_pair->second.neighbor_id;
  parsed_configuration.direction  = periodic_pair->second.direction;

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
    cfd_boundary_conditions.periodic_boundaries.at(boundary_id).neighbor_id;
  parsed_configuration.direction =
    cfd_boundary_conditions.periodic_boundaries.at(boundary_id).direction;

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

  AssertThrow(dem_configuration.direction < dim &&
                cfd_configuration.direction < dim,
              dealii::ExcMessage(
                "Sharp-IB periodic direction must be smaller than the spatial "
                "dimension."));

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

  const double length = domain_upper - domain_lower;
  AssertThrow(std::isfinite(coordinate) && std::isfinite(length) && length > 0.,
              dealii::ExcMessage(
                "Cannot wrap a non-finite point or invalid period."));
  // Keep both boundary faces in place; normalize arbitrarily large excursions.
  if (coordinate < domain_lower)
    wrapped_point[configuration.direction] +=
      std::ceil((domain_lower - coordinate) / length) * length;
  else if (coordinate > domain_upper)
    wrapped_point[configuration.direction] -=
      std::ceil((coordinate - domain_upper) / length) * length;

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
  const double length     = domain_upper - domain_lower;
  AssertThrow(std::isfinite(coordinate) && std::isfinite(support_radius) &&
                support_radius >= 0. && std::isfinite(length) && length > 0.,
              dealii::ExcMessage("Invalid periodic image support or domain."));
  // Enumerate every translated support interval intersecting the domain.
  // Composite constituents may be displaced by more than one period from
  // their reference point, even when the reference point is normalized.
  const double first =
    std::floor((domain_lower - coordinate - support_radius) / length) + 1.;
  const double last =
    std::ceil((domain_upper - coordinate + support_radius) / length) - 1.;
  AssertThrow(
    first > std::numeric_limits<int>::min() &&
      first < std::numeric_limits<int>::max() &&
      last > std::numeric_limits<int>::min() &&
      last < std::numeric_limits<int>::max(),
    dealii::ExcMessage(
      "Periodic image support exceeds the representable image range."));
  for (int image = std::max(1, static_cast<int>(first));
       image <= static_cast<int>(last);
       ++image)
    image_points.push_back(point +
                           static_cast<double>(image) * periodic_offset);
  for (int image = std::min(-1, static_cast<int>(last));
       image >= static_cast<int>(first);
       --image)
    image_points.push_back(point +
                           static_cast<double>(image) * periodic_offset);

  return image_points;
}

template <int dim>
void
SharpIBPeriodicBoundaries<dim>::translate_history(
  std::vector<dealii::Point<dim>> &history,
  const dealii::Tensor<1, dim>    &shift)
{
  for (auto &position : history)
    position += shift;
}

template <int dim>
double
SharpIBPeriodicBoundaries<dim>::enclosing_radius(const Shape<dim> &shape)
{
  return shape.bounding_box_center.norm() +
         shape.bounding_box_half_length.norm();
}

template <int dim>
bool
SharpIBPeriodicBoundaries<dim>::is_periodic_enabled() const
{
  return configuration.enabled;
}

template <int dim>
typename SharpIBPeriodicBoundaries<dim>::PeriodicConfiguration const &
SharpIBPeriodicBoundaries<dim>::get_configuration() const
{
  return configuration;
}

template <int dim>
dealii::Tensor<1, dim>
SharpIBPeriodicBoundaries<dim>::get_periodic_offset() const
{
  return periodic_offset;
}

template <int dim>
std::pair<dealii::types::boundary_id, dealii::types::boundary_id>
SharpIBPeriodicBoundaries<dim>::get_periodic_boundary_ids() const
{
  return {configuration.boundary_0, configuration.boundary_1};
}

template <int dim>
unsigned int
SharpIBPeriodicBoundaries<dim>::get_periodic_direction() const
{
  return configuration.direction;
}

template class SharpIBPeriodicBoundaries<2>;
template class SharpIBPeriodicBoundaries<3>;
