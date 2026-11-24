// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_pipe_and_waveguide_h
#define lethe_grid_pipe_and_waveguide_h

#include <deal.II/base/point.h>
#include <deal.II/base/types.h>

#include <deal.II/grid/tria.h>

#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Generates a 3D mesh of a waveguide embedded in a cylindrical cavity in which a liquid flows
 *
 * The geometry consists of three sections joined along the z-axis:
 * - A bottom cylinder
 * - A central cylinder into which the rectangular-based waveguide is embedded
 * - A top cylinder
 *
 * The cross section of the cylinder is a structured O-grid, designed so that
 * its cells are close to squares:
 * - an inner rectangle around the axis, divided into columns and rows. A band
 *   of rows of width rectangle_width continues the waveguide through the
 *   cylinder, and a strip of rows lies on each side of the band;
 * - a ring of layers of cells between the inner rectangle and the lateral wall.
 *   On the left and right sides, the rows of the band end on the lateral wall
 * at the same height, so that the edges of the band end at the corners
 *   (rectangle_length, 0) and (rectangle_length, rectangle_width) of the
 *   waveguide. The corners of the inner rectangle are joined to the lateral
 * wall at the angle corner_angle from the x-axis, and they are pulled inward so
 *   that the cells of the ring at the corners are less skewed.
 * In the middle section, the waveguide is divided into columns and into the
 * rows of the band, which it shares at the junction with the cylinder. Its
 * number of columns is chosen so that its cells are close to squares.
 *
 * Boundary IDs are assigned as follows:
 * -0: Waveguide port (end face of the waveguide at x = 0)
 * -1: Bottom surface of the cylinder (z = -bottom_height)
 * -2: Top surface of the cylinder (z = middle_height + top_height)
 * -3: Lateral wall of the cylinder
 * -4: Lateral walls of the waveguide
 *
 * Material IDs are assigned as follows:
 * -0: Cylinder
 * -1: Waveguide
 *
 * A CylindricalManifold (manifold id 1) is attached to the faces of the cells
 * of the cylinder that lie on its lateral wall, so that mesh refinement
 * preserves the shape of the cylinder. In the middle section, this includes the
 * face shared by the cylinder and the waveguide, which therefore follows the
 * lateral wall of the cylinder. The other faces and the cells keep the flat
 * manifold, so that the inner rectangle and the walls of the waveguide stay
 * straight, and the interior of the cells next to the lateral wall is
 * interpolated from their faces when they are refined.
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
 */

template <int dim, int spacedim>
class GridPipeAndWaveguide

{
public:
  /**
   * @brief Constructor that parses grid parameters from a colon-separated
   * string, and checks that they define a valid mesh.
   *
   * @param[in] grid_arguments A string with the following format:
   * @code "rectangle_width : rectangle_length : outer_radius :
   * inner_half_width : inner_half_height : corner_angle : corner_pull :
   * n_columns : n_strip_rows : n_band_rows : n_ring_layers : bottom_height :
   * middle_height : top_height : bottom_resolution : middle_resolution :
   * top_resolution" @endcode
   *
   * | #   | Field             | Format        | Description                                                                                         |
   * |-----|-------------------|---------------|-----------------------------------------------------------------------------------------------------|
   * |  0  | rectangle_width   | double        | Waveguide width (along y)                                                                           |
   * |  1  | rectangle_length  | double        | Waveguide length (along x)                                                                          |
   * |  2  | outer_radius      | double        | Radius of the cylinder                                                                              |
   * |  3  | inner_half_width  | double        | Half-width of the inner rectangle (along x)                                                         |
   * |  4  | inner_half_height | double        | Half-height of the inner rectangle (along y), larger than rectangle_width / 2                       |
   * |  5  | corner_angle      | double        | Angle (degrees) from the x-axis of the point of the lateral wall joined to a corner of the rectangle |
   * |  6  | corner_pull       | double        | Fraction of the sides of the corner cells by which the corners of the rectangle are pulled inward   |
   * |  7  | n_columns         | unsigned int  | Number of columns of the inner rectangle                                                            |
   * |  8  | n_strip_rows      | unsigned int  | Number of rows of each strip, on each side of the band                                              |
   * |  9  | n_band_rows       | unsigned int  | Number of rows of the band and of the waveguide                                                     |
   * |  10 | n_ring_layers     | unsigned int  | Number of layers of the ring                                                                        |
   * |  11 | bottom_height     | double        | Height of the bottom cylinder                                                                       |
   * |  12 | middle_height     | double        | Height of the waveguide                                                                             |
   * |  13 | top_height        | double        | Height of the top cylinder                                                                          |
   * |  14 | bottom_resolution | unsigned int  | Number of slices of the extrusion of the bottom cylinder (at least 2)                               |
   * |  15 | middle_resolution | unsigned int  | Number of slices of the extrusion of the middle cylinder (at least 2)                               |
   * |  16 | top_resolution    | unsigned int  | Number of slices of the extrusion of the top cylinder (at least 2)                                  |
   * Example: @code "0.04318 : 0.1 : 0.069 : 0.0363 : 0.0363 : 39.2 : 0.27 : 5 :
   * 1 : 3 : 2 : 0.04932 : 0.08636 : 0.04932 : 5 : 7 : 5" @endcode
   */
  GridPipeAndWaveguide(const std::string &grid_arguments);

  /**
   * @brief Generate the coarse mesh of the cylinder in which the waveguide is
   * embedded.
   *
   * The coarse 2D cross section of the cylinder, and in the middle section of
   * the cylinder and the waveguide, are built by make_cross_section. The three
   * cross sections are extruded along z, and the three sections are merged. The
   * boundary ids and the cylindrical manifold are then set on the merged mesh,
   * which is returned without refinement.
   *
   * @param[out] triangulation The empty triangulation to fill with the mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

  /// Boundary id of the waveguide port (end face of the waveguide at x = 0)
  static constexpr types::boundary_id waveguide_port_boundary_id = 0;

  /// Boundary id of the bottom surface of the cylinder
  static constexpr types::boundary_id bottom_boundary_id = 1;

  /// Boundary id of the top surface of the cylinder
  static constexpr types::boundary_id top_boundary_id = 2;

  /// Boundary id of the lateral wall of the cylinder
  static constexpr types::boundary_id cylinder_wall_boundary_id = 3;

  /// Boundary id of the lateral walls of the waveguide
  static constexpr types::boundary_id waveguide_wall_boundary_id = 4;

  /// Material id of the cells of the cylinder
  static constexpr types::material_id cylinder_material_id = 0;

  /// Material id of the cells of the waveguide
  static constexpr types::material_id waveguide_material_id = 1;

  /// Manifold id of the CylindricalManifold of the cylinder
  static constexpr types::manifold_id cylindrical_manifold_id = 1;

private:
  /**
   * @brief Return the heights of the rows of the inner rectangle, relative to
   * the axis of the cylinder.
   *
   * The n_strip_rows rows of each strip divide it evenly, and so do the
   * n_band_rows rows of the band, so that the edges of the band, at the heights
   * -rectangle_width / 2 and rectangle_width / 2, are row heights.
   *
   * @return The heights of the 2 n_strip_rows + n_band_rows + 1 lines between
   * the rows, from -inner_half_height to inner_half_height.
   */
  std::vector<double>
  row_heights() const;

  /**
   * @brief Return the angle from the x-axis of the point of the lateral wall to
   * which the end of a row line on the right side of the inner rectangle is
   * joined. The left side is the mirror image of the right side.
   *
   * In the band, the line ends on the lateral wall at the same height, so that
   * the cells of the ring in the band are close to rectangles and the edges of
   * the band end at the corners of the waveguide. In the strips, the angle
   * varies linearly with the height, from the angle of the edges of the band to
   * corner_angle at the corners of the rectangle.
   *
   * @param[in] height Height of the row line, relative to the axis.
   *
   * @return The angle in radians.
   */
  double
  wall_angle_of_right_side(const double height) const;

  /**
   * @brief Generate the coarse 2D mesh of the cross section of the cylinder,
   * and of the waveguide if requested.
   *
   * The inner rectangle is divided into n_columns columns and into the rows
   * given by row_heights, and its corners are pulled inward by corner_pull
   * times the sides of their cell. Its boundary is then traversed
   * counterclockwise, and each of its vertices is joined to a point of the
   * lateral wall by n_ring_layers cells of the ring. The waveguide is divided
   * into n_waveguide_columns columns and n_band_rows rows, and its vertices at
   * the junction are the vertices of the band on the lateral wall, so that the
   * two parts share them.
   *
   * The cells of the cylinder and of the waveguide receive their material ids.
   * The mesh has no manifold id: the manifold is set by make_grid on the final
   * 3D mesh.
   *
   * @param[in] center Position of the axis of the cylinder.
   * @param[in] with_waveguide Whether the cells of the waveguide are added.
   * @param[out] triangulation The empty triangulation to fill with the mesh.
   */
  void
  make_cross_section(const Point<2>      &center,
                     const bool           with_waveguide,
                     Triangulation<2, 2> &triangulation) const;

  /// Width of the rectangular waveguide (along y)
  double rectangle_width;

  /// Length of the rectangular waveguide (along x)
  double rectangle_length;

  /// Radius of the cylinder
  double outer_radius;

  /// Half-width of the inner rectangle of the cross section (along x)
  double inner_half_width;

  /// Half-height of the inner rectangle of the cross section (along y)
  double inner_half_height;

  /// Angle (degrees) from the x-axis of the point of the lateral wall joined to
  /// the corner of the inner rectangle
  double corner_angle;

  /// Fraction of the sides of the corner cells of the inner rectangle by which
  /// its corners are pulled inward
  double corner_pull;

  /// Number of columns of the inner rectangle
  unsigned int n_columns;

  /// Number of rows of each strip of the inner rectangle
  unsigned int n_strip_rows;

  /// Number of rows of the band and of the waveguide
  unsigned int n_band_rows;

  /// Number of layers of the ring between the inner rectangle and the wall
  unsigned int n_ring_layers;

  /// Number of columns of the waveguide, chosen so that its cells are close to
  /// squares
  unsigned int n_waveguide_columns;

  /// Height of the bottom cylinder
  double bottom_height;

  /// Height of the middle section (waveguide embedded in the cylinder)
  double middle_height;

  /// Height of the top cylinder
  double top_height;

  /// Number of slices of the extrusion of the bottom cylinder
  unsigned int bottom_resolution;

  /// Number of slices of the extrusion of the middle cylinder
  unsigned int middle_resolution;

  /// Number of slices of the extrusion of the top cylinder
  unsigned int top_resolution;
};


#endif
