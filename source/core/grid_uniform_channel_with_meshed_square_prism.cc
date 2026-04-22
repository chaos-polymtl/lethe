// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_uniform_channel_with_meshed_square_prism.h>
#include <core/utilities.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

template <int dim, int spacedim>
GridUniformChannelWithMeshedSquarePrism<dim, spacedim>::
  GridUniformChannelWithMeshedSquarePrism(const std::string &grid_arguments)
{
  if constexpr (dim == 1 || spacedim == 1)
    {
      AssertThrow(
        false,
        ExcMessage(
          "GridUniformChannelWithMeshedSquarePrism is only supported for <2,2> and <3,3> <dim,spacedim> specializations."));
      return;
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "GridUniformChannelWithMeshedSquarePrism is only supported for <2,2> and <3,3> <dim,spacedim> specializations."));
      return;
    }

  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 5)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory parameters are (bottom left point : top right point : center of the obstacle : inner half-side : outer half-side). The points should be given as x,y and the half-sides should be single numbers. The optional parameters are (rotation_deg : pad_bottom : pad_top : pad_left : pad_right : height : n_slices : colorize )."));
    }

  // Parse bottom_left point
  try
    {
      Tensor<1, 2> bottom_left_coords = value_string_to_tensor<2>(arguments[0]);
      this->bottom_left =
        (dim == 2) ?
          Point<dim>(bottom_left_coords[0], bottom_left_coords[1]) :
          Point<dim>(bottom_left_coords[0], bottom_left_coords[1], 0.0);
    }
  catch (const std::exception &e)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The bottom left point should have 2 components coordinates (x,y format) because the channel is extruded in the z direction."));
    }

  // Parse top_right point
  try
    {
      Tensor<1, 2> top_right_coords = value_string_to_tensor<2>(arguments[1]);
      this->top_right =
        (dim == 2) ? Point<dim>(top_right_coords[0], top_right_coords[1]) :
                     Point<dim>(top_right_coords[0], top_right_coords[1], 0.0);
    }
  catch (const std::exception &e)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The top right point should have 2 components coordinates (x,y format) because the channel is extruded in the z direction."));
    }


  // Parse center point
  try
    {
      Tensor<1, 2> center_coords = value_string_to_tensor<2>(arguments[2]);
      this->center               = (dim == 2) ?
                                     Point<dim>(center_coords[0], center_coords[1]) :
                                     Point<dim>(center_coords[0], center_coords[1], 0.0);
    }
  catch (const std::exception &e)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The center point should have 2 components coordinates (x,y format) because the channel is extruded in the z direction."));
    }


  // Parse half-sides of inner and outer squares
  this->inner_half_side = Utilities::string_to_double(arguments[3]);
  this->outer_half_side = Utilities::string_to_double(arguments[4]);

  AssertThrow(
    inner_half_side > 0.0,
    ExcMessage(
      " The inner_half_side must be greater than 0 to have an object in the channel."));

  AssertThrow(
    outer_half_side >= 1.5 * inner_half_side - 1e-12 * inner_half_side,
    ExcMessage(
      "The outer half-side must be at least 1.5 times the inner half-side, if not the elements in the transition region will be too distorted."));


  // Check if the outer transition square fits in the channel
  AssertThrow((center[0] - outer_half_side >= bottom_left[0]) &&
                (center[0] + outer_half_side <= top_right[0]) &&
                (center[1] - outer_half_side >= bottom_left[1]) &&
                (center[1] + outer_half_side <= top_right[1]),
              ExcMessage(
                "The outer transition square does not fit in the channel."));

  // Parse optional parameters if provided
  this->inner_rotation_angle =
    (arguments.size() > 5) ? Utilities::string_to_double(arguments[5]) : 0.0;

  AssertThrow(
    (inner_rotation_angle >= 0.0) && (inner_rotation_angle < 90.0),
    ExcMessage("The rotation angle needs to be in the range [0, 90) degrees."));

  inner_rotation_angle = inner_rotation_angle * std::numbers::pi / 180.0;

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
  AssertThrow(!((pad_top == 0) && (std::abs((center[1] + outer_half_side) -
                                            top_right[1]) > 1e-12)),
              ExcMessage("Top padding is required because the transition "
                         "region does not reach the top channel boundary."));
  AssertThrow(!((pad_left == 0) && (std::abs((center[0] - outer_half_side) -
                                             bottom_left[0]) > 1e-12)),
              ExcMessage("Left padding is required because the transition "
                         "region does not reach the left channel boundary."));
  AssertThrow(!((pad_right == 0) && (std::abs((center[0] + outer_half_side) -
                                              top_right[0]) > 1e-12)),
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
GridUniformChannelWithMeshedSquarePrism<dim, spacedim>::
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

  Triangulation<2> obstacle_tria;
  // Start with a circular obstacle mesh, which helps define the edges and
  // vertices links. The circular shape is then transformed to the square
  // obstacle by moving the vertices. The outer ring of the circular mesh is
  // used to define the transition region and is attached to the outer square by
  // vertex projection, the inner vertices are shaped to form the inner square.
  // To ensure that the inner square vertices can be easily identified  even
  // when rotated, the initial circular mesh is generated using twice the
  // inner_half_side as radius and every vertex can be identified using a
  // distance criterion.
  GridGenerator::hyper_ball_balanced(obstacle_tria,
                                     center,
                                     2 * inner_half_side);

  // We apply a transformation to shape inner vertices to form the inner square
  GridTools::transform(
    [&](const Point<2> &p) {
      const double dist = p.distance(center);

      // If the vertex is the one in the center, we keep it fixed.
      if (dist <= tol_inner)
        return p;

      // If the vertex is within sqrt(2)*inner_half_side, it cannot be on the
      // outer ring as we defined it to be 2*inner_half_side, so it is an inner
      // vertex. We move it to the inner square by scaling the vector from the
      // center to the vertex. The second condition checks if the vertex is on
      // the diagonal (within a tolerance) in which case we need to scale it to
      // the corner instead of the face of the inner square.
      if ((dist < std::numbers::sqrt2 * inner_half_side + tol_inner) &&
          (std::abs(std::abs(p[0] - center[0]) - std::abs(p[1] - center[1])) <
           tol_inner))
        return center +
               std::numbers::sqrt2 * inner_half_side * (p - center) / dist;
      // This is for non corner inner vertices, we scale them to the face of the
      // inner square.
      if ((dist < std::numbers::sqrt2 * inner_half_side + tol_inner) &&
          (std::abs(std::abs(p[0] - center[0]) - std::abs(p[1] - center[1])) >=
           tol_inner))
        return center + (inner_half_side / dist) * (p - center);
      return p;
    },
    obstacle_tria);

  // We apply the rotation to the whole mesh if a rotation angle is specified.
  if (inner_rotation_angle != 0.0)
    {
      const double cos_a = std::cos(inner_rotation_angle);
      const double sin_a = std::sin(inner_rotation_angle);
      GridTools::transform(
        [&](const Point<2> &p) {
          const double x_rel = p[0] - center[0];
          const double y_rel = p[1] - center[1];
          return Point<2>(center[0] + cos_a * x_rel - sin_a * y_rel,
                          center[1] + sin_a * x_rel + cos_a * y_rel);
        },
        obstacle_tria);
    }

  // Here we pre-compute the 8 projected positions of vertices that we want for
  // the outer square. The outer ring of hyper_ball_balanced has vertices at
  // angles k*pi/4. After rotation they sit at inner_rotation_angle + k*pi/4.
  // The projection maps each direction onto the outer square along the same
  // line. To avoid having voids in the corners we choose the subset of the 4
  // target points that are closest to  corners and change the target projection
  // to attach them there instead.
  std::array<Point<2>, 8> outer_targets;
  for (unsigned int k = 0; k < 8; ++k)
    {
      double theta = inner_rotation_angle + k * std::numbers::pi / 4.0;

      // If the rotation angle is between pi/18 and 4*pi/9, we want to project
      // the vertices at the square middle points to the corners.
      if (inner_rotation_angle >= std::numbers::pi / 18.0 - tol_inner &&
          inner_rotation_angle < 4 * std::numbers::pi / 9.0 - tol_inner)
        {
          // The even k are associated to the vertices at the square middle
          // points.
          if (k % 2 == 0)
            {
              theta = (k + 1) * std::numbers::pi / 4.0;
            }
        }
      // If the rotation angle is between 0 and pi/18 or between 4*pi/9 and
      // pi/2, we want to project the vertices at the square corners to the
      // corners.
      else if (inner_rotation_angle >= 0.0 - tol_inner &&
               inner_rotation_angle < std::numbers::pi / 18.0 - tol_inner)
        {
          // The odd k are associated to the vertices at the square corners.
          if (k % 2 == 1)
            {
              theta = k * std::numbers::pi / 4.0;
            }
        }
      else if (inner_rotation_angle >= 4 * std::numbers::pi / 9.0 - tol_inner &&
               inner_rotation_angle < std::numbers::pi / 2.0 - tol_inner)
        {
          if (k % 2 == 1)
            {
              theta = (k + 2) * std::numbers::pi / 4.0;
            }
        }

      const double cx = std::cos(theta);
      const double cy = std::sin(theta);
      const double scale =
        outer_half_side / std::max(std::abs(cx), std::abs(cy));
      outer_targets[k] = center + Point<2>(scale * cx, scale * cy);
    }

  // Attach each outer ring vertex to its analytical target (closest by angle).
  // Using angle difference avoids the Euclidean distance ambiguity and ensures
  // exact vertex matching with the padding meshes.
  GridTools::transform(
    [&](const Point<2> &p) {
      const double dist = p.distance(center);

      // Check if the vertex is an inside vertex. The outer vertex still lies on
      // the circle defined by the hyper_ball_balanced, so we can use the same
      // distance criterion as before to identify them.
      if (dist <= std::numbers::sqrt2 * inner_half_side + tol_inner)
        return p;

      double   min_angle = -1.0;
      Point<2> closest   = outer_targets[0];
      // We compute the unitary vector from the center to the vertex and do the
      // dot product with the unitary vector from the center to each target. The
      // one with the highest cosine is the closest by angle and we project to
      // that one (when the orientation is superposed, the cosine is 1, when it
      // is orthogonal, the cosine is 0, and when it is opposite, the cosine is
      // -1).
      const Tensor<1, 2> p_normalize = (p - center) / (p - center).norm();
      for (const Point<2> &candidate : outer_targets)
        {
          const Tensor<1, 2> target_normalize =
            (candidate - center) / (candidate - center).norm();
          const double cos_angle = p_normalize * target_normalize;

          if (cos_angle > min_angle)
            {
              min_angle = cos_angle;
              closest   = candidate;
            }
        }
      return closest;
    },
    obstacle_tria);

  // Collect the projected x/y coordinates per outer-square face.
  // Vertices landing exactly on a corner (|cos θ| == |sin θ|) are skipped
  // because they already coincide with the corner of the padding rectangle.
  // This is used to get the padding meshes to have vertices at the exact same
  // positions so that merge_triangulations can merge them without needing to
  // introduce new vertices.
  std::vector<double> top_x, bottom_x, left_y, right_y;
  for (unsigned int k = 0; k < 8; ++k)
    {
      const Point<2> &q      = outer_targets[k];
      const double    theta  = std::atan2(q[1] - center[1], q[0] - center[0]);
      const double    abs_cx = std::abs(std::cos(theta));
      const double    abs_cy = std::abs(std::sin(theta));

      // If the vertex lands on the corner, we do not need to add a breakpoint
      // for the padding meshes because they already have a vertex there.
      if (std::abs(abs_cx - abs_cy) < 1e-12)
        continue;

      // If cos θ is larger than sin θ, then the vertex is either on the left or
      // right face.
      if (abs_cx > abs_cy)
        {
          if (std::cos(theta) > 0.0)
            right_y.push_back(q[1]);
          else
            left_y.push_back(q[1]);
        }
      // If sin θ is larger than cos θ, then the vertex is either on the top or
      // bottom face.
      else
        {
          if (std::sin(theta) > 0.0)
            top_x.push_back(q[0]);
          else
            bottom_x.push_back(q[0]);
        }
    }
  std::ranges::sort(top_x);
  std::ranges::sort(bottom_x);
  std::ranges::sort(left_y);
  std::ranges::sort(right_y);

  // Build a vector of step sizes from sorted breakpoints + outer bounds.
  // This lets subdivided_hyper_rectangle place vertices exactly at the
  // projected positions so merge_triangulations sees coincident nodes.
  // Example:
  //   inner = {x1, x2}, lo = a, hi = b
  //   coords = {a, x1, x2, b}
  //   steps  = {x1-a, x2-x1, b-x2}
  auto make_steps = [](const std::vector<double> &inner_coords,
                       const double               lo,
                       const double               hi) {
    // Create the coordinates list
    std::vector<double> coords = {lo};
    for (const double c : inner_coords)
      coords.push_back(c);
    coords.push_back(hi);

    // Create the steps list
    std::vector<double> steps;
    for (std::size_t i = 1; i < coords.size(); ++i)
      steps.push_back(coords[i] - coords[i - 1]);
    return steps;
  };

  // Similar as above, this lambda function creates a vector of uniform step
  // sizes given the number of steps n and the total length to cover.
  auto uniform = [](const unsigned int n, const double length) {
    return std::vector<double>(n, length / static_cast<double>(n));
  };

  // Bottom pad: columns split at bottom_x projected coords
  Triangulation<2> pad_bottom_tria;
  if (pad_bottom > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_bottom_tria,
      {make_steps(bottom_x,
                  center[0] - outer_half_side,
                  center[0] + outer_half_side),
       uniform(pad_bottom, center[1] - outer_half_side - bottom_left[1])},
      Point<2>(center[0] - outer_half_side, bottom_left[1]),
      Point<2>(center[0] + outer_half_side, center[1] - outer_half_side));

  // Top pad: columns split at top_x projected coords
  Triangulation<2> pad_top_tria;
  if (pad_top > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_top_tria,
      {make_steps(top_x,
                  center[0] - outer_half_side,
                  center[0] + outer_half_side),
       uniform(pad_top, top_right[1] - (center[1] + outer_half_side))},
      Point<2>(center[0] - outer_half_side, center[1] + outer_half_side),
      Point<2>(center[0] + outer_half_side, top_right[1]));

  // Left pad: rows split at left_y projected coords
  Triangulation<2> pad_left_tria;
  if (pad_left > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_left_tria,
      {uniform(pad_left, center[0] - outer_half_side - bottom_left[0]),
       make_steps(left_y,
                  center[1] - outer_half_side,
                  center[1] + outer_half_side)},
      Point<2>(bottom_left[0], center[1] - outer_half_side),
      Point<2>(center[0] - outer_half_side, center[1] + outer_half_side));

  // Right pad: rows split at right_y projected coords
  Triangulation<2> pad_right_tria;
  if (pad_right > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_right_tria,
      {uniform(pad_right, top_right[0] - (center[0] + outer_half_side)),
       make_steps(right_y,
                  center[1] - outer_half_side,
                  center[1] + outer_half_side)},
      Point<2>(center[0] + outer_half_side, center[1] - outer_half_side),
      Point<2>(top_right[0], center[1] + outer_half_side));

  // Corner pads never touch the obstacle tria boundary, so no special split
  // needed
  Triangulation<2> pad_bottom_left_corner_tria;
  if (pad_bottom > 0 && pad_left > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_bottom_left_corner_tria,
      {pad_left, pad_bottom},
      bottom_left,
      Point<2>(center[0] - outer_half_side, center[1] - outer_half_side));

  Triangulation<2> pad_bottom_right_corner_tria;
  if (pad_bottom > 0 && pad_right > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_bottom_right_corner_tria,
      {pad_right, pad_bottom},
      Point<2>(center[0] + outer_half_side, bottom_left[1]),
      Point<2>(top_right[0], center[1] - outer_half_side));

  Triangulation<2> pad_top_left_corner_tria;
  if (pad_top > 0 && pad_left > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_top_left_corner_tria,
      {pad_left, pad_top},
      Point<2>(bottom_left[0], center[1] + outer_half_side),
      Point<2>(center[0] - outer_half_side, top_right[1]));

  Triangulation<2> pad_top_right_corner_tria;
  if (pad_top > 0 && pad_right > 0)
    GridGenerator::subdivided_hyper_rectangle(
      pad_top_right_corner_tria,
      {pad_right, pad_top},
      Point<2>(center[0] + outer_half_side, center[1] + outer_half_side),
      top_right);

  // Merge only non-empty triangulations
  std::vector<const Triangulation<2> *> trias = {&obstacle_tria};
  if (pad_bottom > 0)
    trias.push_back(&pad_bottom_tria);
  if (pad_top > 0)
    trias.push_back(&pad_top_tria);
  if (pad_left > 0)
    trias.push_back(&pad_left_tria);
  if (pad_right > 0)
    trias.push_back(&pad_right_tria);
  if (pad_bottom > 0 && pad_left > 0)
    trias.push_back(&pad_bottom_left_corner_tria);
  if (pad_bottom > 0 && pad_right > 0)
    trias.push_back(&pad_bottom_right_corner_tria);
  if (pad_top > 0 && pad_left > 0)
    trias.push_back(&pad_top_left_corner_tria);
  if (pad_top > 0 && pad_right > 0)
    trias.push_back(&pad_top_right_corner_tria);
  GridGenerator::merge_triangulations(trias, triangulation);

  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids(0);

  // Here we mark obstacle cells with the obstacle material ID. This is done in
  // the unrotated frame by rotating the cell center back to the axis-aligned
  // frame temporarily and checking if it is within the inner square dimensions.
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      const Point<2> c = cell->center();

      // Translate to obstacle-centered coordinates
      const double dx = c[0] - center[0];
      const double dy = c[1] - center[1];

      // Rotate point back to axis-aligned frame
      const double cos_a = std::cos(inner_rotation_angle);
      const double sin_a = std::sin(inner_rotation_angle);

      const double x_local = cos_a * dx + sin_a * dy;
      const double y_local = -sin_a * dx + cos_a * dy;

      // Axis-aligned square test
      const bool inside_obstacle =
        (std::abs(x_local) <= inner_half_side + tol_inner) &&
        (std::abs(y_local) <= inner_half_side + tol_inner);

      cell->set_material_id(inside_obstacle ? obstacle_id : fluid_id);
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
          const Point<2> fc = face->center();
          if (std::abs(fc[0] - bottom_left[0]) < tol_x)
            face->set_boundary_id(0);
          else if (std::abs(fc[0] - top_right[0]) < tol_x)
            face->set_boundary_id(1);
          else if (std::abs(fc[1] - bottom_left[1]) < tol_y)
            face->set_boundary_id(2);
          else if (std::abs(fc[1] - top_right[1]) < tol_y)
            face->set_boundary_id(3);
        }
    }
}

template <>
void
GridUniformChannelWithMeshedSquarePrism<2, 2>::make_grid(
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

  // Attach manifold for id 0 required by deal.II
  TransfiniteInterpolationManifold<2> tfi_manifold;
  tfi_manifold.initialize(triangulation);
  triangulation.set_manifold(0, tfi_manifold);
}


template <>
void
GridUniformChannelWithMeshedSquarePrism<3, 3>::make_grid(
  Triangulation<3, 3> &triangulation)
{
  // Generate the 2D cross-section (geometry + material IDs + boundary IDs)
  Triangulation<2> tria_2D;
  generate_2d_channel_mesh(tria_2D,
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

  // Attach manifold id 0 on the 2D mesh before extrusion.
  // extrude_triangulation queries get_manifold() on the input 2D triangulation.
  TransfiniteInterpolationManifold<2> tfi_manifold_2d;
  tfi_manifold_2d.initialize(tria_2D);
  tria_2D.set_manifold(0, tfi_manifold_2d);

  // Extrude the 2D cross-section along the z-axis. Manifold  and material IDs
  // from the 2D mesh are copied to the lateral faces of the 3D mesh.
  GridGenerator::extrude_triangulation(
    tria_2D, n_slices, height, triangulation, true);

  // Attach manifold for id 0 required by deal.II
  TransfiniteInterpolationManifold<3> tfi_manifold;
  tfi_manifold.initialize(triangulation);
  triangulation.set_manifold(0, tfi_manifold);

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
GridUniformChannelWithMeshedSquarePrism<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridUniformChannelWithMeshedSquarePrism is only supported for <2,2> and <3,3> <dim,spacedim> specializations."));
}

// Explicit template instantiations
template class GridUniformChannelWithMeshedSquarePrism<2, 2>;
template class GridUniformChannelWithMeshedSquarePrism<2, 3>;
template class GridUniformChannelWithMeshedSquarePrism<3, 3>;
