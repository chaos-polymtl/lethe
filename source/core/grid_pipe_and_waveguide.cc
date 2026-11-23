// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_pipe_and_waveguide.h>

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace dealii;

/*
This file creates the mesh of a rectangular waveguide that enters a cylinder
laterally, perpendicular to its axis.

The file is divided into two parts:

1) The function detailed_circle, inspired by GridGenerator::hyper_ball of
deal.II, creates the coarse 2D mesh of the cross section of the cylinder, which
is adapted to be joined to the rectangle of the waveguide on its side.

2) The function make_grid assembles the 3D mesh of the cylinder and of the
waveguide.
*/

// See the Doxygen documentation in the header
template <int dim, int spacedim>
void
GridPipeAndWaveguide<dim, spacedim>::detailed_circle(
  Triangulation<2, 2> &triangulation,
  const Point<2>      &center) const
{
  /*
  The band of width rectangle_width that crosses the circle along the x-axis is
  bounded by the lines y = center[1] +/- half_width.
  - The vertices of the outer circle on its diagonals are at
    (+/- outer_diagonal_coordinate, +/- outer_diagonal_coordinate) from the
    center, and the band crosses the outer circle at
    x = +/- outer_band_crossing_coordinate from the center.
  - The vertices of the inner shape on its diagonals are at
    (+/- inner_diagonal_coordinate, +/- inner_diagonal_coordinate) from the
    center, and the band crosses the inner shape at
    x = +/- inner_band_crossing_coordinate from the center. Only these two
    coordinates depend on the inner shape.
  */
  const double half_width                = rectangle_width / 2.;
  const double outer_diagonal_coordinate = outer_radius / std::numbers::sqrt2;
  const double outer_band_crossing_coordinate =
    std::sqrt(outer_radius * outer_radius - half_width * half_width);
  const double inner_diagonal_coordinate =
    inner_circle ? inner_length / std::numbers::sqrt2 : inner_length;
  const double inner_band_crossing_coordinate =
    inner_circle ?
      std::sqrt(inner_length * inner_length - half_width * half_width) :
      inner_length;

  const std::vector<Point<2>> vertices = {
    // 0-3: Vertices of the outer circle on its diagonals
    center + Point<2>(-outer_diagonal_coordinate, -outer_diagonal_coordinate),
    center + Point<2>(+outer_diagonal_coordinate, -outer_diagonal_coordinate),
    center + Point<2>(-outer_diagonal_coordinate, +outer_diagonal_coordinate),
    center + Point<2>(+outer_diagonal_coordinate, +outer_diagonal_coordinate),
    // 4-7: Vertices of the inner shape on its diagonals
    center + Point<2>(-inner_diagonal_coordinate, +inner_diagonal_coordinate),
    center + Point<2>(+inner_diagonal_coordinate, +inner_diagonal_coordinate),
    center + Point<2>(-inner_diagonal_coordinate, -inner_diagonal_coordinate),
    center + Point<2>(+inner_diagonal_coordinate, -inner_diagonal_coordinate),
    // 8-9: Corners of the rectangle of the waveguide where the band crosses
    // the outer circle on the side of the waveguide
    Point<2>(rectangle_length, 0.),
    Point<2>(rectangle_length, rectangle_width),
    // 10-11: Vertices of the outer circle on the line x = center[0]
    center + Point<2>(0., -outer_radius),
    center + Point<2>(0., +outer_radius),
    // 12-13: Points where the band crosses the outer circle on the side
    // opposite to the waveguide
    center + Point<2>(+outer_band_crossing_coordinate, -half_width),
    center + Point<2>(+outer_band_crossing_coordinate, +half_width),
    // 14-17: Points where the band crosses the inner shape
    center + Point<2>(-inner_band_crossing_coordinate, +half_width),
    center + Point<2>(-inner_band_crossing_coordinate, -half_width),
    center + Point<2>(+inner_band_crossing_coordinate, -half_width),
    center + Point<2>(+inner_band_crossing_coordinate, +half_width),
    // 18-19: Vertices of the inner shape on the line x = center[0]
    center + Point<2>(0., -inner_length),
    center + Point<2>(0., +inner_length),
    // 20-21: Points where the band crosses the line x = center[0]
    center + Point<2>(0., +half_width),
    center + Point<2>(0., -half_width)};

  std::vector<CellData<2>> cells(16, CellData<2>());
  static constexpr int     circle_cell_vertices[16][4] = {
    // We start with the 6 cells inside the inner shape
    {6, 18, 15, 21},
    {15, 21, 14, 20},
    {14, 20, 4, 19},
    {18, 7, 21, 16},
    {21, 16, 20, 17},
    {20, 17, 19, 5},
    // We fill in the outer cells, starting from the top left
    {9, 14, 2, 4},
    {8, 15, 9, 14},
    {0, 6, 8, 15},
    {0, 10, 6, 18},
    {10, 1, 18, 7},
    {7, 1, 16, 12},
    {16, 12, 17, 13},
    {17, 13, 5, 3},
    {19, 5, 11, 3},
    {4, 19, 2, 11}};

  for (unsigned int i = 0; i < 16; ++i)
    {
      for (unsigned int j = 0; j < 4; ++j)
        cells[i].vertices[j] = circle_cell_vertices[i][j];
    }

  triangulation.create_triangulation(vertices,
                                     cells,
                                     SubCellData()); // no boundary information
}



template <int dim, int spacedim>
GridPipeAndWaveguide<dim, spacedim>::GridPipeAndWaveguide(
  const std::string &grid_arguments)
{
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  AssertThrow(
    arguments.size() == 10 || arguments.size() == 11,
    ExcMessage(
      "The pipe_and_waveguide grid needs 10 or 11 grid arguments, but " +
      Utilities::to_string(arguments.size()) +
      " were given. The grid arguments are rectangle_width : "
      "rectangle_length : outer_radius : inner_length : bottom_height : "
      "middle_height : top_height : bottom_resolution : middle_resolution : "
      "top_resolution : inner_circle, where inner_circle is optional and is "
      "true by default."));

  // Lengths
  rectangle_width  = Utilities::string_to_double(arguments[0]);
  rectangle_length = Utilities::string_to_double(arguments[1]);
  outer_radius     = Utilities::string_to_double(arguments[2]);
  inner_length     = Utilities::string_to_double(arguments[3]);
  bottom_height    = Utilities::string_to_double(arguments[4]);
  middle_height    = Utilities::string_to_double(arguments[5]);
  top_height       = Utilities::string_to_double(arguments[6]);

  AssertThrow(
    rectangle_width > 0. && rectangle_length > 0. && outer_radius > 0. &&
      inner_length > 0. && bottom_height > 0. && middle_height > 0. &&
      top_height > 0.,
    ExcMessage(
      "The lengths of the pipe_and_waveguide grid must be positive, but "
      "rectangle_width = " +
      Utilities::to_string(rectangle_width) +
      ", rectangle_length = " + Utilities::to_string(rectangle_length) +
      ", outer_radius = " + Utilities::to_string(outer_radius) +
      ", inner_length = " + Utilities::to_string(inner_length) +
      ", bottom_height = " + Utilities::to_string(bottom_height) +
      ", middle_height = " + Utilities::to_string(middle_height) +
      " and top_height = " + Utilities::to_string(top_height) +
      ". Give a positive value to each of them."));

  // Resolutions. They are read as signed integers so that a negative value is
  // reported instead of being converted to a large unsigned integer.
  const int n_bottom_slices = Utilities::string_to_int(arguments[7]);
  const int n_middle_slices = Utilities::string_to_int(arguments[8]);
  const int n_top_slices    = Utilities::string_to_int(arguments[9]);

  AssertThrow(
    n_bottom_slices >= 2 && n_middle_slices >= 2 && n_top_slices >= 2,
    ExcMessage(
      "bottom_resolution, middle_resolution and top_resolution are the "
      "numbers of slices of the extrusion of the three sections of the "
      "pipe_and_waveguide grid, so each of them must be at least 2, but they "
      "are " +
      Utilities::to_string(n_bottom_slices) + ", " +
      Utilities::to_string(n_middle_slices) + " and " +
      Utilities::to_string(n_top_slices) + ". Set each of them to 2 or more."));

  bottom_resolution = n_bottom_slices;
  middle_resolution = n_middle_slices;
  top_resolution    = n_top_slices;

  // Inner shape. inner_circle is optional and is true by default.
  inner_circle = true;
  if (arguments.size() == 11)
    {
      AssertThrow(arguments[10] == "true" || arguments[10] == "false",
                  ExcMessage(
                    "The last grid argument of the pipe_and_waveguide grid, "
                    "inner_circle, must be true or false, but it is \"" +
                    arguments[10] + "\"."));
      inner_circle = (arguments[10] == "true");
    }

  /*
  Geometric conditions that prevent the generation of inverted cells.

  A band of width rectangle_width crosses the cylinder along the x-axis, in
  line with the waveguide. The inner shape must contain this band and must be
  strictly inside the outer circle.
  - Inner circle of radius inner_length: the points where the edges of the band
    (y = +/- rectangle_width / 2) cross the inner circle must lie between the
    vertices of the inner circle at +/- 45 degrees from the x-axis, hence
    rectangle_width / 2 < inner_length / sqrt(2). Together with
    inner_length < outer_radius, this also ensures that the band crosses the
    outer circle between its vertices at +/- 45 degrees.
  - Inner square of vertices (+/- inner_length, +/- inner_length): the band
    must be narrower than the square, hence rectangle_width / 2 < inner_length,
    and the corners of the square must be inside the outer circle, hence
    inner_length < outer_radius / sqrt(2).
  */
  if (inner_circle)
    {
      AssertThrow(
        inner_length < outer_radius,
        ExcMessage(
          "With an inner circle, inner_length (the radius of the inner "
          "circle) must be smaller than outer_radius, but inner_length = " +
          Utilities::to_string(inner_length) +
          " and outer_radius = " + Utilities::to_string(outer_radius) +
          ". Decrease inner_length or increase outer_radius."));

      AssertThrow(
        rectangle_width < std::numbers::sqrt2 * inner_length,
        ExcMessage(
          "With an inner circle, the inner circle must contain the band of "
          "width rectangle_width that crosses the cylinder, which requires "
          "rectangle_width < sqrt(2) * inner_length, but rectangle_width = " +
          Utilities::to_string(rectangle_width) +
          " and inner_length = " + Utilities::to_string(inner_length) +
          ". Increase inner_length above rectangle_width / sqrt(2) = " +
          Utilities::to_string(rectangle_width / std::numbers::sqrt2) +
          " or decrease rectangle_width."));
    }
  else
    {
      AssertThrow(
        rectangle_width / 2. < inner_length,
        ExcMessage(
          "With an inner square, the inner square must contain the band of "
          "width rectangle_width that crosses the cylinder, which requires "
          "rectangle_width / 2 < inner_length, but rectangle_width = " +
          Utilities::to_string(rectangle_width) +
          " and inner_length = " + Utilities::to_string(inner_length) +
          ". Increase inner_length above rectangle_width / 2 = " +
          Utilities::to_string(rectangle_width / 2.) +
          " or decrease rectangle_width."));

      AssertThrow(
        inner_length < outer_radius / std::numbers::sqrt2,
        ExcMessage(
          "With an inner square, the corners of the inner square must be "
          "inside the outer circle, which requires inner_length < "
          "outer_radius / sqrt(2), but inner_length = " +
          Utilities::to_string(inner_length) +
          " and outer_radius = " + Utilities::to_string(outer_radius) +
          ". Decrease inner_length below outer_radius / sqrt(2) = " +
          Utilities::to_string(outer_radius / std::numbers::sqrt2) +
          " or increase outer_radius."));
    }
}

template <int dim, int spacedim>
void
GridPipeAndWaveguide<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> &triangulation)
{
  // The grid arguments are checked in the constructor

  // Dimension Test
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(false,
                  ExcMessage("GridPipeAndWaveguide is only supported in 3D."));
      return;
    }

  else
    {
      // COARSE 2D CROSS SECTIONS

      // Rectangle of the waveguide, at 0 <= x <= rectangle_length and
      // 0 <= y <= rectangle_width. The waveguide only exists in the middle
      // section, which is why the mesh is built in three sections.
      Triangulation<2> waveguide_rectangle;
      GridGenerator::hyper_rectangle(waveguide_rectangle,
                                     Point<2>(0.0, 0.0),
                                     Point<2>(rectangle_length,
                                              rectangle_width));

      // Cross section of the cylinder, with two vertices at the corners of the
      // rectangle at x = rectangle_length. It is the cross section of the
      // bottom and top sections.
      const Point<2>   center(rectangle_length +
                              std::sqrt(outer_radius * outer_radius -
                                        rectangle_width * rectangle_width / 4),
                            rectangle_width / 2);
      Triangulation<2> disk;
      detailed_circle(disk, center);

      // Cross section of the middle section. The manifold ids are not copied
      // since they are set on the final 3D mesh.
      Triangulation<2> merged;
      GridGenerator::merge_triangulations(
        waveguide_rectangle, disk, merged, 1e-12, false);

      // EXTRUDE

      Triangulation<3> extruded_base;
      Triangulation<3> extruded_middle;
      Triangulation<3> extruded_top;
      GridGenerator::extrude_triangulation(disk,
                                           bottom_resolution,
                                           bottom_height,
                                           extruded_base);
      GridGenerator::extrude_triangulation(merged,
                                           middle_resolution,
                                           middle_height,
                                           extruded_middle);
      GridGenerator::extrude_triangulation(disk,
                                           top_resolution,
                                           top_height,
                                           extruded_top);

      // MERGING THE 3D MESHES

      // Move the cylinders in the right positions to prepare the merge
      Tensor<1, 3> shift_down;
      shift_down[0] = 0.0;
      shift_down[1] = 0.0;
      shift_down[2] = -bottom_height; // The base_cylinder must have its base at
                                      // -bottom_height
      GridTools::shift(shift_down, extruded_base);


      Tensor<1, 3> shift_up;
      shift_up[0] = 0.0;
      shift_up[1] = 0.0;
      shift_up[2] = +middle_height; // The top_cylinder  must have its base at
                                    // +middle_height
      GridTools::shift(shift_up, extruded_top);


      // Merge the 3 parts of the mesh
      std::vector<const Triangulation<3, 3> *> parts = {&extruded_base,
                                                        &extruded_middle,
                                                        &extruded_top};

      // The manifold ids are not copied since they are set on the merged mesh
      GridGenerator::merge_triangulations(parts, triangulation, 1e-12, false);


      // BOUNDARY AND MATERIAL IDs

      /*
      Material ids:
      - cylinder_material_id: cells of the cylinder (x > rectangle_length)
      - waveguide_material_id: cells of the waveguide (x <= rectangle_length)

      Boundary ids:
      - waveguide_port_boundary_id: end face of the waveguide (x = 0)
      - waveguide_wall_boundary_id: other boundary faces of the waveguide
      - bottom_boundary_id: bottom surface of the cylinder
        (z = -bottom_height)
      - top_boundary_id: top surface of the cylinder
        (z = middle_height + top_height)
      - cylinder_wall_boundary_id: lateral wall of the cylinder

      The material id of a cell is set first, since it is used to distinguish
      the walls of the waveguide from the walls of the cylinder.
      */

      for (auto &cell : triangulation.active_cell_iterators())
        {
          const bool cell_is_in_waveguide =
            cell->center()[0] <= rectangle_length;
          cell->set_material_id(cell_is_in_waveguide ? waveguide_material_id :
                                                       cylinder_material_id);

          for (const auto f : cell->face_indices())
            {
              if (!cell->face(f)->at_boundary())
                continue;

              const Point<3> face_center = cell->face(f)->center();
              if (face_center[0] <= 1e-6)
                cell->face(f)->set_boundary_id(waveguide_port_boundary_id);
              else if (cell_is_in_waveguide)
                cell->face(f)->set_boundary_id(waveguide_wall_boundary_id);
              else if (face_center[2] <= -bottom_height + 1e-6)
                cell->face(f)->set_boundary_id(bottom_boundary_id);
              else if (face_center[2] >= middle_height + top_height - 1e-6)
                cell->face(f)->set_boundary_id(top_boundary_id);
              else
                cell->face(f)->set_boundary_id(cylinder_wall_boundary_id);
            }
        }


      // MANIFOLD

      /*
      The cylindrical manifold is attached to the cells of the cylinder that
      have a face on its lateral wall, that is a face whose vertices are all at
      a distance outer_radius from the axis of the cylinder. This includes the
      cells of the middle section that are next to the waveguide: their face
      shared with the waveguide is an interior face of the mesh, but it lies on
      the lateral wall of the cylinder. The manifold ids are set on the merged
      mesh so that they are the same in the three sections. The cells of the
      waveguide keep the flat manifold, so that the walls of the waveguide stay
      planar.
      - Inner circle: the manifold id is set on the whole cell, so that its
        edges on the inner circle, whose vertices are all at a distance
        inner_length from the axis, are also arcs of a circle.
      - Inner square: the manifold id is only set on the face on the lateral
        wall. The edges of the cell on the inner square join vertices at
        different distances from the axis, and the cylindrical manifold would
        bend them, whereas the sides of the square must stay straight. The
        interior of the cell is interpolated from its faces when it is refined.
      */
      const double distance_to_axis_tolerance = 1e-8 * outer_radius;
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          if (cell->material_id() != cylinder_material_id)
            continue;

          for (const auto f : cell->face_indices())
            {
              bool face_is_on_lateral_wall = true;
              for (const auto v : cell->face(f)->vertex_indices())
                {
                  const Point<3> &vertex = cell->face(f)->vertex(v);
                  const double    distance_to_axis =
                    center.distance(Point<2>(vertex[0], vertex[1]));
                  if (std::abs(distance_to_axis - outer_radius) >
                      distance_to_axis_tolerance)
                    face_is_on_lateral_wall = false;
                }

              if (face_is_on_lateral_wall)
                {
                  if (inner_circle)
                    cell->set_all_manifold_ids(cylindrical_manifold_id);
                  else
                    cell->face(f)->set_all_manifold_ids(
                      cylindrical_manifold_id);
                  break;
                }
            }
        }

      const Point<3>     center_3d(center[0], center[1], 0.0);
      const Tensor<1, 3> axis({0.0, 0.0, 1.0});
      triangulation.set_manifold(cylindrical_manifold_id,
                                 CylindricalManifold<3>(axis, center_3d));
    }
}

// Explicit instantiations
template class GridPipeAndWaveguide<2, 2>;
template class GridPipeAndWaveguide<2, 3>;
template class GridPipeAndWaveguide<3, 3>;
