// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_pipe_and_waveguide.h>

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

using namespace dealii;

/*
This file creates the mesh of a rectangular waveguide that enters a cylinder
laterally, perpendicular to its axis.

The function make_cross_section creates the coarse 2D mesh of the cross section:
a structured O-grid of the cylinder around an inner rectangle, and in the middle
section the waveguide. The function make_grid extrudes the cross sections of
the three sections and assembles the 3D mesh.
*/

namespace
{
  /**
   * @brief Return whether a quadrilateral cell is valid, that is whether the
   * Jacobian of its bilinear map is positive at its four vertices.
   *
   * @param[in] vertices The vertices of the mesh.
   * @param[in] cell The cell, with its vertices in the lexicographic order of
   * deal.II: (0, 0), (1, 0), (0, 1) and (1, 1) in the reference cell.
   *
   * @return True if the Jacobian is positive at the four vertices.
   */
  bool
  cell_is_valid(const std::vector<Point<2>> &vertices, const CellData<2> &cell)
  {
    // Edges of the cell along the two directions of the reference cell
    const Tensor<1, 2> bottom_edge =
      vertices[cell.vertices[1]] - vertices[cell.vertices[0]];
    const Tensor<1, 2> top_edge =
      vertices[cell.vertices[3]] - vertices[cell.vertices[2]];
    const Tensor<1, 2> left_edge =
      vertices[cell.vertices[2]] - vertices[cell.vertices[0]];
    const Tensor<1, 2> right_edge =
      vertices[cell.vertices[3]] - vertices[cell.vertices[1]];

    // At each vertex, the Jacobian is the cross product of the two edges that
    // start from it
    const std::vector<std::pair<Tensor<1, 2>, Tensor<1, 2>>> vertex_edges = {
      {bottom_edge, left_edge},
      {bottom_edge, right_edge},
      {top_edge, left_edge},
      {top_edge, right_edge}};
    for (const auto &[first_edge, second_edge] : vertex_edges)
      if (first_edge[0] * second_edge[1] - first_edge[1] * second_edge[0] <= 0.)
        return false;
    return true;
  }
} // namespace

template <int dim, int spacedim>
GridPipeAndWaveguide<dim, spacedim>::GridPipeAndWaveguide(
  const std::string &grid_arguments)
{
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  AssertThrow(
    arguments.size() == 17,
    ExcMessage(
      "The pipe_and_waveguide grid needs 17 grid arguments, but " +
      Utilities::to_string(arguments.size()) +
      " were given. The grid arguments are rectangle_width : "
      "rectangle_length : outer_radius : inner_half_width : "
      "inner_half_height : corner_angle : corner_pull : n_columns : "
      "n_strip_rows : n_band_rows : n_ring_layers : bottom_height : "
      "middle_height : top_height : bottom_resolution : middle_resolution : "
      "top_resolution."));

  // Lengths and parameters of the inner rectangle
  rectangle_width   = Utilities::string_to_double(arguments[0]);
  rectangle_length  = Utilities::string_to_double(arguments[1]);
  outer_radius      = Utilities::string_to_double(arguments[2]);
  inner_half_width  = Utilities::string_to_double(arguments[3]);
  inner_half_height = Utilities::string_to_double(arguments[4]);
  corner_angle      = Utilities::string_to_double(arguments[5]);
  corner_pull       = Utilities::string_to_double(arguments[6]);
  bottom_height     = Utilities::string_to_double(arguments[11]);
  middle_height     = Utilities::string_to_double(arguments[12]);
  top_height        = Utilities::string_to_double(arguments[13]);

  AssertThrow(
    rectangle_width > 0. && rectangle_length > 0. && outer_radius > 0. &&
      inner_half_width > 0. && inner_half_height > 0. && bottom_height > 0. &&
      middle_height > 0. && top_height > 0.,
    ExcMessage(
      "The lengths of the pipe_and_waveguide grid must be positive, but "
      "rectangle_width = " +
      Utilities::to_string(rectangle_width) +
      ", rectangle_length = " + Utilities::to_string(rectangle_length) +
      ", outer_radius = " + Utilities::to_string(outer_radius) +
      ", inner_half_width = " + Utilities::to_string(inner_half_width) +
      ", inner_half_height = " + Utilities::to_string(inner_half_height) +
      ", bottom_height = " + Utilities::to_string(bottom_height) +
      ", middle_height = " + Utilities::to_string(middle_height) +
      " and top_height = " + Utilities::to_string(top_height) +
      ". Give a positive value to each of them."));

  // Numbers of cells and of slices. They are read as signed integers so that
  // a negative value is reported instead of being converted to a large
  // unsigned integer.
  const int n_columns_argument     = Utilities::string_to_int(arguments[7]);
  const int n_strip_rows_argument  = Utilities::string_to_int(arguments[8]);
  const int n_band_rows_argument   = Utilities::string_to_int(arguments[9]);
  const int n_ring_layers_argument = Utilities::string_to_int(arguments[10]);
  const int n_bottom_slices        = Utilities::string_to_int(arguments[14]);
  const int n_middle_slices        = Utilities::string_to_int(arguments[15]);
  const int n_top_slices           = Utilities::string_to_int(arguments[16]);

  AssertThrow(
    n_columns_argument >= 1 && n_strip_rows_argument >= 1 &&
      n_band_rows_argument >= 1 && n_ring_layers_argument >= 1,
    ExcMessage(
      "n_columns, n_strip_rows, n_band_rows and n_ring_layers are numbers of "
      "cells of the cross section of the pipe_and_waveguide grid, so each of "
      "them must be at least 1, but they are " +
      Utilities::to_string(n_columns_argument) + ", " +
      Utilities::to_string(n_strip_rows_argument) + ", " +
      Utilities::to_string(n_band_rows_argument) + " and " +
      Utilities::to_string(n_ring_layers_argument) +
      ". Set each of them to 1 or more."));

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

  n_columns         = n_columns_argument;
  n_strip_rows      = n_strip_rows_argument;
  n_band_rows       = n_band_rows_argument;
  n_ring_layers     = n_ring_layers_argument;
  bottom_resolution = n_bottom_slices;
  middle_resolution = n_middle_slices;
  top_resolution    = n_top_slices;

  /*
  Geometric conditions on the inner rectangle:
  - it contains the band of width rectangle_width that continues the waveguide,
    with a strip on each side of the band;
  - its corners are inside the cylinder;
  - its corners are joined to the lateral wall between the edges of the band
    and the y-axis, so that the cells of the ring follow each other around the
    rectangle.
  The cells of the cross section are also checked, in make_cross_section.
  */
  AssertThrow(
    inner_half_height > rectangle_width / 2.,
    ExcMessage(
      "The inner rectangle of the pipe_and_waveguide grid must contain the "
      "band of width rectangle_width that continues the waveguide, which "
      "requires inner_half_height > rectangle_width / 2, but "
      "inner_half_height = " +
      Utilities::to_string(inner_half_height) +
      " and rectangle_width = " + Utilities::to_string(rectangle_width) +
      ". Increase inner_half_height above " +
      Utilities::to_string(rectangle_width / 2.) + "."));

  AssertThrow(
    inner_half_width * inner_half_width +
        inner_half_height * inner_half_height <
      outer_radius * outer_radius,
    ExcMessage(
      "The corners of the inner rectangle of the pipe_and_waveguide grid must "
      "be inside the cylinder, which requires sqrt(inner_half_width^2 + "
      "inner_half_height^2) < outer_radius, but inner_half_width = " +
      Utilities::to_string(inner_half_width) +
      ", inner_half_height = " + Utilities::to_string(inner_half_height) +
      " and outer_radius = " + Utilities::to_string(outer_radius) +
      ". Decrease inner_half_width or inner_half_height."));

  const double band_edge_angle =
    std::asin(rectangle_width / 2. / outer_radius) * 180. / std::numbers::pi;
  AssertThrow(
    corner_angle > band_edge_angle && corner_angle < 90.,
    ExcMessage(
      "corner_angle, the angle from the x-axis of the point of the lateral "
      "wall joined to a corner of the inner rectangle of the "
      "pipe_and_waveguide grid, must be between the angle of the edges of "
      "the band, " +
      Utilities::to_string(band_edge_angle) + ", and 90 degrees, but it is " +
      Utilities::to_string(corner_angle) + "."));

  AssertThrow(
    corner_pull >= 0. && corner_pull < 0.5,
    ExcMessage(
      "corner_pull, the fraction of the sides of the corner cells by which the "
      "corners of the inner rectangle of the pipe_and_waveguide grid are "
      "pulled inward, must be in [0, 0.5), but it is " +
      Utilities::to_string(corner_pull) + "."));

  /*
  The number of columns of the waveguide is chosen so that its cells are close
  to squares. At the junction, the lateral wall bulges into the waveguide, so
  the cells of the last column must be longer than the depth of this bulge.
  */
  n_waveguide_columns = static_cast<unsigned int>(
    std::max(1L,
             std::lround(rectangle_length * n_band_rows / rectangle_width)));
  const double wall_bulge_depth =
    outer_radius - std::sqrt(outer_radius * outer_radius -
                             rectangle_width * rectangle_width / 4.);
  AssertThrow(
    rectangle_length / n_waveguide_columns > wall_bulge_depth,
    ExcMessage(
      "The cells of the waveguide of the pipe_and_waveguide grid must be "
      "longer than the depth " +
      Utilities::to_string(wall_bulge_depth) +
      " by which the lateral wall bulges into the waveguide, but they are " +
      Utilities::to_string(rectangle_length / n_waveguide_columns) +
      " long. Decrease n_band_rows or increase rectangle_length."));
}

template <int dim, int spacedim>
std::vector<double>
GridPipeAndWaveguide<dim, spacedim>::row_heights() const
{
  const double half_band_width = rectangle_width / 2.;
  const double strip_row_height =
    (inner_half_height - half_band_width) / n_strip_rows;
  const double band_row_height = rectangle_width / n_band_rows;

  std::vector<double> heights;
  heights.reserve(2 * n_strip_rows + n_band_rows + 1);
  for (unsigned int s = 0; s < n_strip_rows; ++s)
    heights.push_back(-inner_half_height + s * strip_row_height);
  for (unsigned int b = 0; b < n_band_rows; ++b)
    heights.push_back(-half_band_width + b * band_row_height);
  for (unsigned int s = 0; s < n_strip_rows; ++s)
    heights.push_back(half_band_width + s * strip_row_height);
  heights.push_back(inner_half_height);
  return heights;
}

template <int dim, int spacedim>
double
GridPipeAndWaveguide<dim, spacedim>::wall_angle_of_right_side(
  const double height) const
{
  const double half_band_width = rectangle_width / 2.;
  const double tolerance       = 1e-12 * outer_radius;

  // In the band, the row line ends on the lateral wall at the same height
  if (std::abs(height) <= half_band_width + tolerance)
    return std::asin(std::clamp(height / outer_radius, -1., 1.));

  // In the strips, the angle varies linearly with the height, from the edge of
  // the band to the corner of the rectangle
  const double band_edge_angle = std::asin(half_band_width / outer_radius);
  const double corner_angle_in_radians = corner_angle * std::numbers::pi / 180.;
  const double fraction                = (std::abs(height) - half_band_width) /
                          (inner_half_height - half_band_width);
  const double angle =
    band_edge_angle + fraction * (corner_angle_in_radians - band_edge_angle);
  return height > 0. ? angle : -angle;
}

template <int dim, int spacedim>
void
GridPipeAndWaveguide<dim, spacedim>::make_cross_section(
  const Point<2>      &center,
  const bool           with_waveguide,
  Triangulation<2, 2> &triangulation) const
{
  const std::vector<double> heights    = row_heights();
  const unsigned int        n_rows     = heights.size() - 1;
  const double corner_angle_in_radians = corner_angle * std::numbers::pi / 180.;
  const unsigned int n_boundary_vertices = 2 * (n_columns + n_rows);

  std::vector<Point<2>>    vertices;
  std::vector<CellData<2>> cells;
  vertices.reserve((n_columns + 1) * (n_rows + 1) +
                   n_ring_layers * n_boundary_vertices +
                   n_waveguide_columns * (n_band_rows + 1));
  cells.reserve(n_columns * n_rows + n_ring_layers * n_boundary_vertices +
                n_waveguide_columns * n_band_rows);

  // INNER RECTANGLE

  // The vertex (i, j) is at the intersection of the column line i and of the
  // row line j, and its index is j (n_columns + 1) + i
  for (unsigned int j = 0; j <= n_rows; ++j)
    for (unsigned int i = 0; i <= n_columns; ++i)
      vertices.push_back(
        center +
        Point<2>(-inner_half_width +
                   2. * inner_half_width * static_cast<double>(i) / n_columns,
                 heights[j]));

  // The four corners are pulled inward along the diagonal of their cell, by
  // the fraction corner_pull of the sides of the cell. The angle of the
  // rectangle at its corners then exceeds 90 degrees, and the two cells of the
  // ring at each corner share less than 270 degrees, which makes them less
  // skewed.
  const double column_width = 2. * inner_half_width / n_columns;
  const double strip_row_height =
    (inner_half_height - rectangle_width / 2.) / n_strip_rows;
  for (const unsigned int i : {0U, n_columns})
    for (const unsigned int j : {0U, n_rows})
      {
        const double direction_x = (i == 0) ? 1. : -1.;
        const double direction_y = (j == 0) ? 1. : -1.;
        vertices[j * (n_columns + 1) + i] +=
          Tensor<1, 2>({direction_x * corner_pull * column_width,
                        direction_y * corner_pull * strip_row_height});
      }

  for (unsigned int j = 0; j < n_rows; ++j)
    for (unsigned int i = 0; i < n_columns; ++i)
      {
        CellData<2> cell;
        cell.vertices[0] = j * (n_columns + 1) + i;
        cell.vertices[1] = j * (n_columns + 1) + i + 1;
        cell.vertices[2] = (j + 1) * (n_columns + 1) + i;
        cell.vertices[3] = (j + 1) * (n_columns + 1) + i + 1;
        cell.material_id = cylinder_material_id;
        cells.push_back(cell);
      }

  // RING

  // The boundary of the inner rectangle is traversed counterclockwise from its
  // bottom left corner: bottom side from left to right, right side from bottom
  // to top, top side from right to left and left side from top to bottom. Each
  // of its vertices is joined to the point of the lateral wall at the angle
  // wall_angles from the x-axis.
  std::vector<unsigned int> boundary_vertices;
  std::vector<double>       wall_angles;
  boundary_vertices.reserve(n_boundary_vertices);
  wall_angles.reserve(n_boundary_vertices);
  const double top_and_bottom_wall_angle =
    std::numbers::pi - 2. * corner_angle_in_radians;
  for (unsigned int i = 0; i < n_columns; ++i)
    {
      boundary_vertices.push_back(i);
      wall_angles.push_back(std::numbers::pi + corner_angle_in_radians +
                            top_and_bottom_wall_angle * static_cast<double>(i) /
                              n_columns);
    }
  for (unsigned int j = 0; j < n_rows; ++j)
    {
      boundary_vertices.push_back(j * (n_columns + 1) + n_columns);
      wall_angles.push_back(wall_angle_of_right_side(heights[j]));
    }
  for (unsigned int i = n_columns; i > 0; --i)
    {
      boundary_vertices.push_back(n_rows * (n_columns + 1) + i);
      wall_angles.push_back(corner_angle_in_radians +
                            top_and_bottom_wall_angle *
                              static_cast<double>(n_columns - i) / n_columns);
    }
  for (unsigned int j = n_rows; j > 0; --j)
    {
      boundary_vertices.push_back(j * (n_columns + 1));
      wall_angles.push_back(std::numbers::pi -
                            wall_angle_of_right_side(heights[j]));
    }

  // The vertex k of the layer l of the ring lies on the segment between the
  // boundary vertex k and its point of the lateral wall, at the fraction
  // l / n_ring_layers of the segment. The vertices of the layer 0 are the
  // boundary vertices, and those of the last layer are on the lateral wall.
  std::vector<std::vector<unsigned int>> ring_vertices(n_ring_layers + 1,
                                                       boundary_vertices);
  for (unsigned int l = 1; l <= n_ring_layers; ++l)
    for (unsigned int k = 0; k < n_boundary_vertices; ++k)
      {
        const Point<2> boundary_point = vertices[boundary_vertices[k]];
        const Point<2> wall_point =
          center + Point<2>(outer_radius * std::cos(wall_angles[k]),
                            outer_radius * std::sin(wall_angles[k]));
        const double fraction = static_cast<double>(l) / n_ring_layers;
        ring_vertices[l][k]   = vertices.size();
        vertices.push_back(l == n_ring_layers ?
                             wall_point :
                             boundary_point +
                               fraction * (wall_point - boundary_point));
      }

  // The cell k of the layer l joins the vertices k and k + 1 of the layers l
  // and l + 1. Its first reference direction goes from the vertex k + 1 to the
  // vertex k, and its second one goes outward, so that it is positively
  // oriented.
  for (unsigned int l = 0; l < n_ring_layers; ++l)
    for (unsigned int k = 0; k < n_boundary_vertices; ++k)
      {
        const unsigned int next = (k + 1) % n_boundary_vertices;
        CellData<2>        cell;
        cell.vertices[0] = ring_vertices[l][next];
        cell.vertices[1] = ring_vertices[l][k];
        cell.vertices[2] = ring_vertices[l + 1][next];
        cell.vertices[3] = ring_vertices[l + 1][k];
        cell.material_id = cylinder_material_id;
        cells.push_back(cell);
      }

  // WAVEGUIDE

  // The vertex (i, j) of the waveguide is at the intersection of the column
  // line i and of the row line j of the band. The column lines are evenly
  // spaced between the port (x = 0) and rectangle_length, except the last
  // one: its vertices are the vertices of the band on the lateral wall, which
  // the waveguide shares with the ring.
  if (with_waveguide)
    {
      std::vector<std::vector<unsigned int>> waveguide_vertices(
        n_band_rows + 1, std::vector<unsigned int>(n_waveguide_columns + 1));
      for (unsigned int j = 0; j <= n_band_rows; ++j)
        {
          // The row line j of the band is the row line n_strip_rows + j of the
          // rectangle. On the left side, which is traversed last, the boundary
          // vertex of the row line j_row has the index
          // 2 n_columns + n_rows + (n_rows - j_row).
          const unsigned int row_line = n_strip_rows + j;
          const unsigned int wall_vertex =
            ring_vertices[n_ring_layers]
                         [2 * n_columns + n_rows + (n_rows - row_line)];
          waveguide_vertices[j][n_waveguide_columns] = wall_vertex;

          const double row_y = vertices[wall_vertex][1];
          for (unsigned int i = 0; i < n_waveguide_columns; ++i)
            {
              waveguide_vertices[j][i] = vertices.size();
              vertices.emplace_back(rectangle_length * static_cast<double>(i) /
                                      n_waveguide_columns,
                                    row_y);
            }
        }

      for (unsigned int j = 0; j < n_band_rows; ++j)
        for (unsigned int i = 0; i < n_waveguide_columns; ++i)
          {
            CellData<2> cell;
            cell.vertices[0] = waveguide_vertices[j][i];
            cell.vertices[1] = waveguide_vertices[j][i + 1];
            cell.vertices[2] = waveguide_vertices[j + 1][i];
            cell.vertices[3] = waveguide_vertices[j + 1][i + 1];
            cell.material_id = waveguide_material_id;
            cells.push_back(cell);
          }
    }

  for (const auto &cell : cells)
    AssertThrow(
      cell_is_valid(vertices, cell),
      ExcMessage(
        "The grid arguments of the pipe_and_waveguide grid give an inverted "
        "cell in the cross section, with inner_half_width = " +
        Utilities::to_string(inner_half_width) +
        ", inner_half_height = " + Utilities::to_string(inner_half_height) +
        ", corner_angle = " + Utilities::to_string(corner_angle) +
        " and corner_pull = " + Utilities::to_string(corner_pull) +
        ". Move the corners of the inner rectangle closer to the direction "
        "of corner_angle, or decrease corner_pull."));

  triangulation.create_triangulation(vertices, cells, SubCellData());
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

      // The axis of the cylinder is placed so that the lateral wall passes
      // through the corners (rectangle_length, 0) and
      // (rectangle_length, rectangle_width) of the waveguide. The waveguide
      // only exists in the middle section, which is why the mesh is built in
      // three sections.
      const Point<2>   center(rectangle_length +
                              std::sqrt(outer_radius * outer_radius -
                                        rectangle_width * rectangle_width / 4.),
                            rectangle_width / 2.);
      Triangulation<2> cylinder_cross_section;
      Triangulation<2> middle_cross_section;
      make_cross_section(center, false, cylinder_cross_section);
      make_cross_section(center, true, middle_cross_section);

      // EXTRUDE

      Triangulation<3> extruded_base;
      Triangulation<3> extruded_middle;
      Triangulation<3> extruded_top;
      GridGenerator::extrude_triangulation(cylinder_cross_section,
                                           bottom_resolution,
                                           bottom_height,
                                           extruded_base);
      GridGenerator::extrude_triangulation(middle_cross_section,
                                           middle_resolution,
                                           middle_height,
                                           extruded_middle);
      GridGenerator::extrude_triangulation(cylinder_cross_section,
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

      // The manifold ids are not copied since they are set on the merged mesh.
      // The material ids, set on the cross sections, are kept by the extrusion
      // and the merge.
      GridGenerator::merge_triangulations(parts, triangulation, 1e-12, false);


      // BOUNDARY IDs

      /*
      - waveguide_port_boundary_id: end face of the waveguide (x = 0)
      - waveguide_wall_boundary_id: other boundary faces of the waveguide
      - bottom_boundary_id: bottom surface of the cylinder
        (z = -bottom_height)
      - top_boundary_id: top surface of the cylinder
        (z = middle_height + top_height)
      - cylinder_wall_boundary_id: lateral wall of the cylinder
      */
      for (const auto &cell : triangulation.active_cell_iterators())
        {
          const bool cell_is_in_waveguide =
            cell->material_id() == waveguide_material_id;

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
      The cylindrical manifold is attached to the faces of the cells of the
      cylinder that lie on its lateral wall, that is the faces whose vertices
      are all at a distance outer_radius from the axis of the cylinder. This
      includes the faces shared by the cylinder and the waveguide in the middle
      section: they are interior faces of the mesh, but they lie on the lateral
      wall of the cylinder. The manifold ids are set on the merged mesh so that
      they are the same in the three sections.

      Only the faces receive the manifold id. The other edges of these cells
      join vertices at different distances from the axis, and the cylindrical
      manifold would bend them, whereas the inner rectangle, the layers of the
      ring and the walls of the waveguide must stay straight. The interior of
      the cells is interpolated from their faces when they are refined.
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
                  cell->face(f)->set_all_manifold_ids(cylindrical_manifold_id);
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
