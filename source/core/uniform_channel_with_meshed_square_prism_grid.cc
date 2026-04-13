// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/uniform_channel_with_meshed_square_prism_grid.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

template <int dim, int spacedim>
UniformChannelWithMeshedSquarePrismGrid<dim, spacedim>::
  UniformChannelWithMeshedSquarePrismGrid(
    const std::string &grid_arguments)
{
  if constexpr (dim == 1 || spacedim == 1)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with rotated meshed rectangle obstacle is only supported in 2D and 3D space."));
      return;
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with rotated meshed rectangle obstacle is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 5)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory parameters are (bottom left point : top right point : center of the obstacle : inner half-side : outer half-side). The points should be given as x,y and the half-sides should be single numbers. The optional parameters are (rotation_deg : pad bottom : pad top : pad left : pad right : height : n_slices : colorize )."));
    }

  // Parse bottom_left point
  std::stringstream   bottom_left_stream(arguments[0]);
  std::vector<double> bottom_left_coords;
  std::string         coord_str;
  while (getline(bottom_left_stream, coord_str, ','))
    {
      bottom_left_coords.push_back(Utilities::string_to_double(coord_str));
    }

  if (bottom_left_coords.size() != 2)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The bottom left point should have 2 coordinates (x,y format) because the channel is extruded in the z direction."));
    }

  this->bottom_left =
    (dim == 2) ? Point<dim>(bottom_left_coords[0], bottom_left_coords[1]) :
                 Point<dim>(bottom_left_coords[0], bottom_left_coords[1], 0.0);

  // Parse top_right point
  std::stringstream   top_right_stream(arguments[1]);
  std::vector<double> top_right_coords;
  while (getline(top_right_stream, coord_str, ','))
    {
      top_right_coords.push_back(Utilities::string_to_double(coord_str));
    }

  if (top_right_coords.size() != 2)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The top right point should have 2 coordinates (x,y format) because the channel is extruded in the z direction."));
    }

  this->top_right = (dim == 2) ?
                      Point<dim>(top_right_coords[0], top_right_coords[1]) :
                      Point<dim>(top_right_coords[0], top_right_coords[1], 0.0);

  // Parse center point
  std::stringstream   center_stream(arguments[2]);
  std::vector<double> center_coords;
  while (getline(center_stream, coord_str, ','))
    {
      center_coords.push_back(Utilities::string_to_double(coord_str));
    }

  if (center_coords.size() != 2)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The center point should have 2 coordinates (x,y format) because the channel is extruded in the z direction."));
    }

  this->center = (dim == 2) ?
                   Point<dim>(center_coords[0], center_coords[1]) :
                   Point<dim>(center_coords[0], center_coords[1], 0.0);

  // Parse half-sides of inner and outer squares
  this->inner_half_side = Utilities::string_to_double(arguments[3]);
  this->outer_half_side = Utilities::string_to_double(arguments[4]);

  AssertThrow(inner_half_side > 0.0 && outer_half_side >= inner_half_side,
              ExcMessage("Require outer_half_side >= inner_half_side > 0."));

  // Check if the outer transition square fits in the channel
  AssertThrow((center[0] - outer_half_side >= bottom_left[0]) &&
                (center[0] + outer_half_side <= top_right[0]) &&
                (center[1] - outer_half_side >= bottom_left[1]) &&
                (center[1] + outer_half_side <= top_right[1]),
              ExcMessage(
                "The outer transition square does not fit in the channel."));

  // Parse optional parameters if provided
  this->inner_rotation_angle =
    (arguments.size() > 5) ? std::remainder(Utilities::string_to_double(arguments[5]), 90.0) *
                               numbers::PI / 180.0 :
                             0.0;
  this->pad_bottom =
    (arguments.size() > 6) ? Utilities::string_to_int(arguments[6]) : 0;
  this->pad_top =
    (arguments.size() > 7) ? Utilities::string_to_int(arguments[7]) : 0;
  this->pad_left =
    (arguments.size() > 8) ? Utilities::string_to_int(arguments[8]) : 0;
  this->pad_right =
    (arguments.size() > 9) ? Utilities::string_to_int(arguments[9]) : 0;

  // Padding is required in each direction where the transition square does not
  // reach the channel boundary (i.e., there is a gap to fill with rectangles).
  AssertThrow(!((pad_bottom == 0) && (std::abs((center[1] - outer_half_side) -
                                               bottom_left[1]) > 1e-12)),
              ExcMessage("Bottom padding is required because the transition "
                         "region does not reach the bottom channel boundary."));
  AssertThrow(!((pad_top == 0) &&
                (std::abs((center[1] + outer_half_side) - top_right[1]) > 1e-12)),
              ExcMessage("Top padding is required because the transition "
                         "region does not reach the top channel boundary."));
  AssertThrow(!((pad_left == 0) && (std::abs((center[0] - outer_half_side) -
                                             bottom_left[0]) > 1e-12)),
              ExcMessage("Left padding is required because the transition "
                         "region does not reach the left channel boundary."));
  AssertThrow(!((pad_right == 0) &&
                (std::abs((center[0] + outer_half_side) - top_right[0]) > 1e-12)),
              ExcMessage("Right padding is required because the transition "
                         "region does not reach the right channel boundary."));

  this->height =
    (arguments.size() > 10) ? Utilities::string_to_double(arguments[10]) : 1.0;
  this->n_slices =
    (arguments.size() > 11) ? Utilities::string_to_int(arguments[11]) : 2;
  this->colorize =
    (arguments.size() > 12 && arguments[12] == "true") ? true : false;
}

template <int dim, int spacedim>
void
UniformChannelWithMeshedSquarePrismGrid<dim, spacedim>::
  generate_2d_channel_mesh(Triangulation<2>  &triangulation,
                           const Point<2>    &bottom_left,
                           const Point<2>    &top_right,
                           const Point<2>    &center,
                           const double       inner_half_side,
                           const double       outer_half_side,
                           const double       inner_rotation_angle,
                           const unsigned int pad_bottom,
                           const unsigned int pad_top,
                           const unsigned int pad_left,
                           const unsigned int pad_right,
                           const bool         colorize)
{
  const types::material_id fluid_id    = 0;
  const types::material_id obstacle_id = 1;
  const double             tol_inner   = 1e-12 * inner_half_side;

  //Verify that the inner square can be created with the desired rotation without exceeding the outer square.
  const double max_diagonal_half =
    inner_half_side * std::numbers::sqrt2 * std::abs(std::sin(inner_rotation_angle));
  AssertThrow(max_diagonal_half <= outer_half_side,
              ExcMessage(
                "The inner square cannot be inscribed in the outer square with the desired rotation. Increase the outer_half_side or decrease the rotation angle."));

  // Create a balanced disk triangulation that the inner square obstacle is
  // inscribed in. This provides a structured inner mesh before square shaping.
  Triangulation<2> obstacle_tria;
  GridGenerator::hyper_ball_balanced(obstacle_tria, center, 2 * inner_half_side);

  // Move inside vertices to form the inner square.
  GridTools::transform(
    [&](const Point<2> &p) {
      const double dist = p.distance(center);
      if (dist <= tol_inner)
        return p;

      // Vertices that lie inside the circle and that are on diagonals needs to be moved by sqrt2 * inner_half_side to reach the square boundary.
      if (( dist < std::numbers::sqrt2 * inner_half_side) &&
          (std::abs(std::abs(p[0] - center[0]) - std::abs(p[1] - center[1])) <
           tol_inner))
        {
          return center +
                 std::numbers::sqrt2 * inner_half_side * (p - center) / dist;
        }
      // Vertices that lie inside the circle and that are not on diagonals needs to be moved by inner_half_side to reach the square boundary.
      if ((dist < std::numbers::sqrt2 * inner_half_side) &&
          (std::abs(std::abs(p[0] - center[0]) - std::abs(p[1] - center[1])) >=
           tol_inner))
        {
          return center + (inner_half_side / dist) * (p - center);
        }

      return p;
    },
    obstacle_tria);

  // Rotate the obstacle triangulation once by the requested angle.
  if (inner_rotation_angle != 0.0)
    {
      const double cos_angle = std::cos(inner_rotation_angle);
      const double sin_angle = std::sin(inner_rotation_angle);
      GridTools::transform(
        [&](const Point<2> &p) {
          const double x_rel = p[0] - center[0];
          const double y_rel = p[1] - center[1];
          return Point<2>(center[0] + cos_angle * x_rel - sin_angle * y_rel,
                          center[1] + sin_angle * x_rel + cos_angle * y_rel);
        },
        obstacle_tria);
    }

  // Scale vertices on the circle to the outer square boundary max length.
  GridTools::transform(
    [&](const Point<2> &p) {
      const double dist = p.distance(center);
      if (std::abs(dist - 2 * inner_half_side) < tol_inner){
        std::cout << "Vertex " << p << " and move to " <<
                  center +  std::numbers::sqrt2 *  (outer_half_side / dist) * (p - center) << std::endl;
        return center + std::numbers::sqrt2 * outer_half_side * (p - center) / dist;
      }
      return p;
    },
    obstacle_tria);

  // Associate each outer-ring vertex to the closest target point on the outer
  // square boundary.
  const std::array<Point<2>, 8> target_points = {{
    center + Point<2>(outer_half_side, outer_half_side),
    center + Point<2>(-outer_half_side, outer_half_side),
    center + Point<2>(outer_half_side, -outer_half_side),
    center + Point<2>(-outer_half_side, -outer_half_side),
    center + Point<2>(outer_half_side, 0.0),
    center + Point<2>(-outer_half_side, 0.0),
    center + Point<2>(0.0, outer_half_side),
    center + Point<2>(0.0, -outer_half_side),
              }};
  GridTools::transform(
    [&](const Point<2> &p) {
      const double dist = p.distance(center);
      if ( dist <= std::numbers::sqrt2 * inner_half_side)
        return p;

      double   min_dist_sq = std::numeric_limits<double>::max();
      Point<2> closest     = target_points[0];
      for (const Point<2> &candidate : target_points)
        {
          const double dx   = p[0] - candidate[0];
          const double dy   = p[1] - candidate[1];
          const double d_sq = dx * dx + dy * dy;
          if (d_sq < min_dist_sq)
            {
              min_dist_sq = d_sq;
              closest     = candidate;
            }
        }
        std::cout << "Vertex " << p << " is moved to " << closest << " with distance "
                  << std::sqrt(min_dist_sq) << std::endl;
      return closest;
    },
    obstacle_tria);


  // Create the padding with up to 8 rectangles that fill the gap between the
  // transition square and the channel boundary. The four sides and four
  // corners are meshed independently.
  Triangulation<2> pad_bottom_tria;
  if (pad_bottom > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_tria,
        {2, pad_bottom},
        Point<2>(center[0] - outer_half_side, bottom_left[1]),
        Point<2>(center[0] + outer_half_side, center[1] - outer_half_side));
    }

  Triangulation<2> pad_top_tria;
  if (pad_top > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_tria,
        {2, pad_top},
        Point<2>(center[0] - outer_half_side, center[1] + outer_half_side),
        Point<2>(center[0] + outer_half_side, top_right[1]));
    }

  Triangulation<2> pad_left_tria;
  if (pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_left_tria,
        {pad_left, 2},
        Point<2>(bottom_left[0], center[1] - outer_half_side),
        Point<2>(center[0] - outer_half_side, center[1] + outer_half_side));
    }

  Triangulation<2> pad_right_tria;
  if (pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_right_tria,
        {pad_right, 2},
        Point<2>(center[0] + outer_half_side, center[1] - outer_half_side),
        Point<2>(top_right[0], center[1] + outer_half_side));
    }

  Triangulation<2> pad_bottom_left_corner_tria;
  if (pad_bottom > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_left_corner_tria,
        {pad_left, pad_bottom},
        bottom_left,
        Point<2>(center[0] - outer_half_side, center[1] - outer_half_side));
    }

  Triangulation<2> pad_bottom_right_corner_tria;
  if (pad_bottom > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_right_corner_tria,
        {pad_right, pad_bottom},
        Point<2>(center[0] + outer_half_side, bottom_left[1]),
        Point<2>(top_right[0], center[1] - outer_half_side));
    }

  Triangulation<2> pad_top_left_corner_tria;
  if (pad_top > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_left_corner_tria,
        {pad_left, pad_top},
        Point<2>(bottom_left[0], center[1] + outer_half_side),
        Point<2>(center[0] - outer_half_side, top_right[1]));
    }

  Triangulation<2> pad_top_right_corner_tria;
  if (pad_top > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_right_corner_tria,
        {pad_right, pad_top},
        Point<2>(center[0] + outer_half_side, center[1] + outer_half_side),
        top_right);
    }

  // Merge all sub-triangulations into the final channel mesh
  GridGenerator::merge_triangulations({&obstacle_tria,
                                       &pad_bottom_tria,
                                       &pad_top_tria,
                                       &pad_left_tria,
                                       &pad_right_tria,
                                       &pad_bottom_left_corner_tria,
                                       &pad_bottom_right_corner_tria,
                                       &pad_top_left_corner_tria,
                                       &pad_top_right_corner_tria},
                                      triangulation);


  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids(0);

  // Mark obstacle cells with material_id=1 based on the inner square region. If we are inside the inner square, we are in the obstacle.
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      const Point<2> cell_center = cell->center();
      if ((std::abs(cell_center[0] - center[0]) <= inner_half_side) &&
          (std::abs(cell_center[1] - center[1]) <= inner_half_side + 1e-12))
        {
          cell->set_material_id(obstacle_id);
        }
      else
        {
          cell->set_material_id(fluid_id);
        }
      }

  // Assign boundary IDs following the subdivided_hyper_rectangle convention:
  //   0: left (-x),  1: right (+x),  2: bottom (-y),  3: top (+y)
  if (colorize)
    {
      const double tol_x = 1e-10 * (top_right[0] - bottom_left[0]);
      const double tol_y = 1e-10 * (top_right[1] - bottom_left[1]);

      for (const auto &face : triangulation.active_face_iterators())
        {
          if (!face->at_boundary())
            continue;

          const Point<2> face_center = face->center();

          if (std::abs(face_center[0] - bottom_left[0]) < tol_x)
            face->set_boundary_id(0);
          else if (std::abs(face_center[0] - top_right[0]) < tol_x)
            face->set_boundary_id(1);
          else if (std::abs(face_center[1] - bottom_left[1]) < tol_y)
            face->set_boundary_id(2);
          else if (std::abs(face_center[1] - top_right[1]) < tol_y)
            face->set_boundary_id(3);
        }
    }
}

template <>
void
UniformChannelWithMeshedSquarePrismGrid<2, 2>::make_grid(
  Triangulation<2, 2> &triangulation)
{
  generate_2d_channel_mesh(triangulation,
                           bottom_left,
                           top_right,
                           center,
                           inner_half_side,
                           outer_half_side,
                           inner_rotation_angle,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);
}


template <>
void
UniformChannelWithMeshedSquarePrismGrid<3, 3>::make_grid(
  Triangulation<3, 3> &triangulation)
{
  // Generate the 2D cross-section (geometry + material IDs + boundary IDs)
  Triangulation<2> tria_D;
  generate_2d_channel_mesh(tria_D,
                           Point<2>(bottom_left[0], bottom_left[1]),
                           Point<2>(top_right[0], top_right[1]),
                           Point<2>(center[0], center[1]),
                           inner_half_side,
                           outer_half_side,
                            inner_rotation_angle,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);

  GridGenerator::extrude_triangulation(
    tria_D, n_slices, height, triangulation, true);

  // Set the boundary ids for the extruded top and bottom faces.
  if (colorize)
    {
      for (const auto &face : triangulation.active_face_iterators())
        {
          if (!face->at_boundary())
            continue;

          const Point<3> face_center = face->center();

          if (std::abs(face_center[2] - bottom_left[2]) < 1e-10 * height)
            face->set_boundary_id(4); // bottom
          else if (std::abs(face_center[2] - (bottom_left[2] + height)) <
                   1e-10 * height)
            face->set_boundary_id(5); // top
        }
    }
}

// Fallback make_grid definition for unsupported template parameters.
template <int dim, int spacedim>
void
UniformChannelWithMeshedSquarePrismGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "UniformChannelWithMeshedSquarePrismGrid is only supported for <2,2> and <3,3> <dim,spacedim> specializations."));
}

// Explicit template instantiations
template class UniformChannelWithMeshedSquarePrismGrid<2, 2>;
template class UniformChannelWithMeshedSquarePrismGrid<2, 3>;
template class UniformChannelWithMeshedSquarePrismGrid<3, 3>;
