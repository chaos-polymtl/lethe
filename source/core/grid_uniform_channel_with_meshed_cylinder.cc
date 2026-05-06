// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_uniform_channel_with_meshed_cylinder.h>
#include <core/utilities.h>

#include <numbers>

template <int dim, int spacedim>
GridUniformChannelWithMeshedCylinder<dim, spacedim>::
  GridUniformChannelWithMeshedCylinder(const std::string &grid_arguments)
{
  if constexpr (dim == 1 || spacedim == 1)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with meshed cylinder mesh is only supported in 2D and 3D space."));
      return;
    }
  else if constexpr (dim == 2 && spacedim == 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The uniform channel with meshed cylinder mesh is only supported in 3D space with 3D elements."));
      return;
    }

  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 5)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory uniform channel with meshed cylinder parameters are (bottom left point : top right point : center of the cylinder : inner radius : outer radius). The points should be given as x,y and the radii should be given as a single number. The optional parameters are (pad bottom : pad top : pad left : pad right : height : n_slices : use_transfinite_region : colorize). The padding parameters should be given as a single integer, the height should be given as a single double, the n_slices should be given as a single integer, use_transfinite_region should be given as true/false, and the colorize parameter should be given as true/false."));
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

  // Parse inner_radius
  this->inner_radius = Utilities::string_to_double(arguments[3]);

  // Parse outer_radius
  this->outer_radius = Utilities::string_to_double(arguments[4]);

  // Check if the cylinder fits in the rectangle
  AssertThrow((center[0] - outer_radius >= bottom_left[0]) &&
                (center[0] + outer_radius <= top_right[0]) &&
                (center[1] - outer_radius >= bottom_left[1]) &&
                (center[1] + outer_radius <= top_right[1]),
              ExcMessage(
                "The cylinder outer radius does not fit in the channel."));

  // Parse optional parameters if provided. If not provided, they will be set to
  // default values.
  this->pad_bottom =
    (arguments.size() > 5) ? Utilities::string_to_int(arguments[5]) : 0;
  this->pad_top =
    (arguments.size() > 6) ? Utilities::string_to_int(arguments[6]) : 0;
  this->pad_left =
    (arguments.size() > 7) ? Utilities::string_to_int(arguments[7]) : 0;
  this->pad_right =
    (arguments.size() > 8) ? Utilities::string_to_int(arguments[8]) : 0;

  // Padding is required in each direction where the transition square does not
  // reach the channel boundary (i.e., there is a gap to fill with rectangles).
  AssertThrow(!((pad_bottom == 0) && (std::abs((center[1] - outer_radius) -
                                               bottom_left[1]) > 1e-12)),
              ExcMessage("Bottom padding is required because the transition "
                         "region does not reach the bottom channel boundary."));
  AssertThrow(!((pad_top == 0) &&
                (std::abs((center[1] + outer_radius) - top_right[1]) > 1e-12)),
              ExcMessage("Top padding is required because the transition "
                         "region does not reach the top channel boundary."));
  AssertThrow(!((pad_left == 0) && (std::abs((center[0] - outer_radius) -
                                             bottom_left[0]) > 1e-12)),
              ExcMessage("Left padding is required because the transition "
                         "region does not reach the left channel boundary."));
  AssertThrow(!((pad_right == 0) &&
                (std::abs((center[0] + outer_radius) - top_right[0]) > 1e-12)),
              ExcMessage("Right padding is required because the transition "
                         "region does not reach the right channel boundary."));

  this->height =
    (arguments.size() > 9) ? Utilities::string_to_double(arguments[9]) : 1.0;
  this->n_slices =
    (arguments.size() > 10) ? Utilities::string_to_int(arguments[10]) : 2.;
  this->use_transfinite_region =
    (arguments.size() > 11 && arguments[11] == "true") ? true : false;
  this->colorize =
    (arguments.size() > 12 && arguments[12] == "true") ? true : false;
}


template <int dim, int spacedim>
void
GridUniformChannelWithMeshedCylinder<dim, spacedim>::generate_2d_channel_mesh(
  Triangulation<2>  &triangulation,
  const Point<2>    &bottom_left,
  const Point<2>    &top_right,
  const Point<2>    &center,
  const double       inner_radius,
  const double       outer_radius,
  const unsigned int pad_bottom,
  const unsigned int pad_top,
  const unsigned int pad_left,
  const unsigned int pad_right,
  const bool         colorize)
{
  const types::material_id fluid_material_id = 0;
  const types::material_id solid_material_id = 1;
  const types::manifold_id tfi_manifold_id   = 0;
  const types::manifold_id polar_manifold_id = 1;

  // Create the inner cylinder (balanced disk)
  Triangulation<2> cylinder_tria;
  GridGenerator::hyper_ball_balanced(cylinder_tria, center, inner_radius);

  // Create the transition region between the cylinder and the rectangular
  // channel. A hyper_shell with 8 segments produces vertices at 45-degree
  // intervals. The diagonal vertices are moved outward by a factor of sqrt(2)
  // so that the outer boundary forms a square of half-side outer_radius.
  Triangulation<2> box_tria;
  GridGenerator::hyper_shell(box_tria, center, inner_radius, outer_radius, 8);
  for (const auto &cell : box_tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_cell; ++v)
        {
          Point<2>    &vertex = cell->vertex(v);
          const double dist   = vertex.distance(center);
          // Check if the vertex is on the outer ring AND on a diagonal
          if ((std::abs(dist - outer_radius) < 1e-12 * outer_radius) &&
              (std::abs(std::abs(vertex[0] - center[0]) -
                        std::abs(vertex[1] - center[1]))) <
                1e-12 * outer_radius)
            {
              vertex = center + std::numbers::sqrt2 * (vertex - center);
            }
        }
    }

  // Create the padding with up to 8 rectangles that fill the gap between the
  // transition square and the channel boundary. The four sides and four
  // corners are meshed independently.
  Triangulation<2> pad_bottom_tria;
  if (pad_bottom > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_tria,
        {2, pad_bottom},
        Point<2>(center[0] - outer_radius, bottom_left[1]),
        Point<2>(center[0] + outer_radius, center[1] - outer_radius));
    }

  Triangulation<2> pad_top_tria;
  if (pad_top > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_tria,
        {2, pad_top},
        Point<2>(center[0] - outer_radius, center[1] + outer_radius),
        Point<2>(center[0] + outer_radius, top_right[1]));
    }

  Triangulation<2> pad_left_tria;
  if (pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_left_tria,
        {pad_left, 2},
        Point<2>(bottom_left[0], center[1] - outer_radius),
        Point<2>(center[0] - outer_radius, center[1] + outer_radius));
    }

  Triangulation<2> pad_right_tria;
  if (pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_right_tria,
        {pad_right, 2},
        Point<2>(center[0] + outer_radius, center[1] - outer_radius),
        Point<2>(top_right[0], center[1] + outer_radius));
    }

  Triangulation<2> pad_bottom_left_corner_tria;
  if (pad_bottom > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_left_corner_tria,
        {pad_left, pad_bottom},
        bottom_left,
        Point<2>(center[0] - outer_radius, center[1] - outer_radius));
    }

  Triangulation<2> pad_bottom_right_corner_tria;
  if (pad_bottom > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_bottom_right_corner_tria,
        {pad_right, pad_bottom},
        Point<2>(center[0] + outer_radius, bottom_left[1]),
        Point<2>(top_right[0], center[1] - outer_radius));
    }

  Triangulation<2> pad_top_left_corner_tria;
  if (pad_top > 0 && pad_left > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_left_corner_tria,
        {pad_left, pad_top},
        Point<2>(bottom_left[0], center[1] + outer_radius),
        Point<2>(center[0] - outer_radius, top_right[1]));
    }

  Triangulation<2> pad_top_right_corner_tria;
  if (pad_top > 0 && pad_right > 0)
    {
      GridGenerator::subdivided_hyper_rectangle(
        pad_top_right_corner_tria,
        {pad_right, pad_top},
        Point<2>(center[0] + outer_radius, center[1] + outer_radius),
        top_right);
    }

  // Merge only non-empty triangulations
  std::vector<const Triangulation<2> *> trias = {&cylinder_tria, &box_tria};
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


  // Assign material and manifold IDs:
  //  - id 0 (TFI) on all cells inside the transition region and on all faces
  //  not on the inner cylinder boundary
  //  - id 1 (polar/cylindrical) on inner-cylinder boundary faces
  triangulation.reset_all_manifolds();
  triangulation.set_all_manifold_ids(tfi_manifold_id);

  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // The inner cylinder is marked with the polar manifold ID for the
      // material and all faces that have all vertices on the inner circle are
      // marked with the polar manifold ID for the manifold.
      if (cell->center().distance(center) < inner_radius)
        {
          cell->set_material_id(solid_material_id);
          for (const auto &face : cell->face_iterators())
            {
              bool all_vertices_on_circle = true;
              for (unsigned int v = 0; v < GeometryInfo<2>::vertices_per_face;
                   ++v)
                {
                  if (std::abs(face->vertex(v).distance(center) -
                               inner_radius) > 1e-10 * inner_radius)
                    {
                      all_vertices_on_circle = false;
                      break;
                    }
                }
              if (all_vertices_on_circle)
                face->set_all_manifold_ids(polar_manifold_id);
            }
        }
      // If outside of the transition region, we can mark cells as flat manifold
      // and fluid material.
      else if (std::abs(cell->center()[0] - center[0]) >
                 outer_radius + 1e-10 * outer_radius ||
               std::abs(cell->center()[1] - center[1]) >
                 outer_radius + 1e-10 * outer_radius)
        {
          cell->set_all_manifold_ids(numbers::flat_manifold_id);
          cell->set_material_id(fluid_material_id);
        }
      // Every other cells stay with TFI manifold and gets the fluid material
      // ID.
      else
        {
          cell->set_material_id(fluid_material_id);
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
GridUniformChannelWithMeshedCylinder<2, 2>::make_grid(
  Triangulation<2, 2> &triangulation)
{
  generate_2d_channel_mesh(triangulation,
                           bottom_left,
                           top_right,
                           center,
                           inner_radius,
                           outer_radius,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);

  // Attach manifold objects for proper refinement behavior
  FlatManifold<2, 2> flat_manifold;
  triangulation.set_manifold(numbers::flat_manifold_id, flat_manifold);

  PolarManifold<2, 2> polar_manifold(center);
  triangulation.set_manifold(1, polar_manifold);

  triangulation.set_manifold(0, FlatManifold<2, 2>());
  if (use_transfinite_region)
    {
      TransfiniteInterpolationManifold<2> tfi_manifold;
      tfi_manifold.initialize(triangulation);
      triangulation.set_manifold(0, tfi_manifold);
    }
}


template <>
void
GridUniformChannelWithMeshedCylinder<3, 3>::make_grid(
  Triangulation<3, 3> &triangulation)
{
  // Generate the 2D cross-section (geometry + manifold IDs + boundary IDs)
  Triangulation<2> tria_2D;
  generate_2d_channel_mesh(tria_2D,
                           Point<2>(bottom_left[0], bottom_left[1]),
                           Point<2>(top_right[0], top_right[1]),
                           Point<2>(center[0], center[1]),
                           inner_radius,
                           outer_radius,
                           pad_bottom,
                           pad_top,
                           pad_left,
                           pad_right,
                           colorize);

  // Extrude the 2D cross-section along the z-axis. Manifold  and material IDs
  // from the 2D mesh are copied to the lateral faces of the 3D mesh.
  GridGenerator::extrude_triangulation(
    tria_2D, n_slices, height, triangulation, true);

  // Attach 3D manifold objects. The manifold IDs (0 for TFI, 1 for
  // cylindrical) were inherited from the 2D mesh during extrusion.
  FlatManifold<3, 3> flat_manifold;
  triangulation.set_manifold(numbers::flat_manifold_id, flat_manifold);

  const Tensor<1, 3>           direction{{0.0, 0.0, 1.0}};
  const CylindricalManifold<3> cylindrical_manifold(direction, center);
  triangulation.set_manifold(1, cylindrical_manifold);

  triangulation.set_manifold(0, FlatManifold<3, 3>());
  if (use_transfinite_region)
    {
      TransfiniteInterpolationManifold<3> tfi_manifold;
      tfi_manifold.initialize(triangulation);
      triangulation.set_manifold(0, tfi_manifold);
    }

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

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
GridUniformChannelWithMeshedCylinder<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridUniformChannelWithMeshedCylinder is only supported for <2,2> and <3,3> <dim,spacedim> specializations."));
}

// Explicit template instantiations
template class GridUniformChannelWithMeshedCylinder<2, 2>;
template class GridUniformChannelWithMeshedCylinder<2, 3>;
template class GridUniformChannelWithMeshedCylinder<3, 3>;
