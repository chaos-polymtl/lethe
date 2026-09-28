// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/immersed_solid_classifier.h>

#include <deal.II/base/exceptions.h>

#include <boost/geometry/index/predicates.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <algorithm>

using namespace dealii;

namespace
{
  /// Relative enlargement of the ball enclosing a cell. The bounding box given
  /// by the mapping is built from its support points, which can slightly
  /// underestimate the extent of a curved cell.
  constexpr double cell_ball_safety_factor = 1.01;
} // namespace

template <int dim>
CellPhase
ImmersedSolidClassifier<dim>::classify_cell(
  const typename Triangulation<dim>::cell_iterator &cell,
  const Mapping<dim>                               &mapping,
  std::vector<unsigned int>                        &solids) const
{
  const BoundingBox<dim> cell_box        = mapping.get_bounding_box(cell);
  const auto &[lower_point, upper_point] = cell_box.get_boundary_points();
  const double radius =
    cell_ball_safety_factor * 0.5 * lower_point.distance(upper_point);

  return classify_ball(cell_box.center(), radius, solids);
}

template <int dim>
ShapesImmersedSolidClassifier<dim>::ShapesImmersedSolidClassifier(
  const std::vector<std::shared_ptr<Shape<dim>>> &solid_shapes)
  : shapes(solid_shapes)
{
  std::vector<std::pair<BoundingBox<dim>, unsigned int>> bounding_boxes;
  for (unsigned int s = 0; s < shapes.size(); ++s)
    {
      const std::shared_ptr<Shape<dim>> &shape = shapes[s];

      AssertThrow(shape != nullptr,
                  ExcMessage("The shape of an immersed solid is not "
                             "defined."));
      AssertThrow(shape->additional_info_on_shape != "iges",
                  ExcMessage("Shapes defined from IGES files are not "
                             "supported since they do not define an inside "
                             "and an outside."));

      // The bounding box of a sphere is built from its exact bounding ball.
      // The signed distance at the center of the sphere is minus its radius,
      // layer thickening included.
      if (dynamic_cast<const Sphere<dim> *>(shape.get()) != nullptr)
        {
          const Point<dim> center = shape->get_position();
          const double     radius = -shape->value(center);
          AssertThrow(radius > 0.,
                      ExcMessage("The radius of a spherical immersed solid, "
                                 "layer thickening included, must be strictly "
                                 "positive."));

          Point<dim> lower_point(center);
          Point<dim> upper_point(center);
          for (unsigned int d = 0; d < dim; ++d)
            {
              lower_point[d] -= radius;
              upper_point[d] += radius;
            }
          bounding_boxes.emplace_back(
            BoundingBox<dim>(std::make_pair(lower_point, upper_point)), s);
        }
      else
        unbounded_solids.push_back(s);
    }

  bounded_solids_tree = pack_rtree(bounding_boxes);
}

template <int dim>
CellPhase
ShapesImmersedSolidClassifier<dim>::classify_ball(
  const Point<dim>          &center,
  const double               radius,
  std::vector<unsigned int> &solids) const
{
  // Gather the candidate solids: the solids whose bounding box intersects the
  // box enclosing the ball and the solids without bounding box.
  solids.clear();
  Point<dim> lower_point(center);
  Point<dim> upper_point(center);
  for (unsigned int d = 0; d < dim; ++d)
    {
      lower_point[d] -= radius;
      upper_point[d] += radius;
    }
  const BoundingBox<dim> ball_box(std::make_pair(lower_point, upper_point));
  bounded_solids_tree.query(
    boost::geometry::index::intersects(ball_box),
    boost::make_function_output_iterator(
      [&solids](const std::pair<BoundingBox<dim>, unsigned int> &leaf) {
        solids.push_back(leaf.second);
      }));
  solids.insert(solids.end(), unbounded_solids.begin(), unbounded_solids.end());

  // Sort the candidates so that, when the ball lies in overlapping solids, the
  // solid with the lowest index is reported, as in the sharp-edge immersed
  // boundary solver.
  std::sort(solids.begin(), solids.end());

  // Since the signed distance is 1-Lipschitz, a solid whose distance to the
  // center is larger than the radius cannot intersect the ball, and a solid
  // whose distance is smaller than minus the radius contains the ball. The
  // candidates that may cut the ball are compacted at the front of the vector.
  unsigned int n_cutting_solids = 0;
  for (unsigned int i = 0; i < solids.size(); ++i)
    {
      const unsigned int solid           = solids[i];
      const double       signed_distance = shapes[solid]->value(center);
      if (signed_distance < -radius)
        {
          solids.assign(1, solid);
          return CellPhase::solid;
        }
      if (signed_distance <= radius)
        solids[n_cutting_solids++] = solid;
    }
  solids.resize(n_cutting_solids);

  return (n_cutting_solids == 0) ? CellPhase::fluid : CellPhase::cut;
}

template <int dim>
void
ShapesImmersedSolidClassifier<dim>::fill_fluid_indicator(
  const std::vector<unsigned int> &cutting_solids,
  const std::vector<Point<dim>>   &points,
  std::vector<double>             &fluid_indicator) const
{
  fluid_indicator.resize(points.size());
  for (unsigned int q = 0; q < points.size(); ++q)
    {
      // A point with a non-positive signed distance is inside a solid, as in
      // the sharp-edge immersed boundary solver.
      double indicator = 1.;
      for (const unsigned int solid : cutting_solids)
        {
          if (shapes[solid]->value(points[q]) <= 0.)
            {
              indicator = 0.;
              break;
            }
        }
      fluid_indicator[q] = indicator;
    }
}

template class ImmersedSolidClassifier<2>;
template class ImmersedSolidClassifier<3>;
template class ShapesImmersedSolidClassifier<2>;
template class ShapesImmersedSolidClassifier<3>;
