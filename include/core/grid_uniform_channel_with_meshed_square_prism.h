// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_uniform_channel_with_meshed_square_prism_h
#define lethe_grid_uniform_channel_with_meshed_square_prism_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <cmath>
#include <sstream>

using namespace dealii;

/**
 * @brief Generates a rectangular channel mesh with a meshed square obstacle
 * in 2D (extruded to a square prism in 3D). By default, this geometry defines
 * two material IDs. The obstacle is assigned the material_id = 1 and the
 * remaining domain is assigned material_id = 0.
 *
 * The mesh consists of three nested regions around the obstacle center:
 * - An inner square of half-side @p inner_half_side (material_id = 1).
 * - A transition ring between @p inner_half_side and @p outer_half_side.
 * - Rectangular padding regions that fill the gap between the transition
 *    square and the channel boundary defined by @p bottom_left / @p top_right.
 *
 * The number of subdivisions in each padding direction (bottom, top, left,
 * right) can be controlled independently. In 3D, the 2D cross-section is
 * extruded along the z-axis with a configurable height and number of slices.
 * The slices refer to the number of layers in the z direction, so the total
 * number of cells in the z direction will be n_slices - 1. It follows that the
 * minimum number of slices is 2, which corresponds to a single layer of cells
 * in the z direction. Note that it is always assumed that the channel bottom is
 * at z=0 and the channel top is at z=height.
 *
 */
template <int dim, int spacedim>
class GridUniformChannelWithMeshedSquarePrism
{
public:
  /**
   * @brief Constructor that parses geometry parameters from a colon-separated
   * string.
   *
   * @param[in] grid_arguments A colon-separated string with the following
   * fields (coordinates are comma-separated):
   *
   * | #  | Field           | Format     | Required | Description                                                          |
   * |----|-----------------|------------|----------|----------------------------------------------------------------------|
   * |  0 | bottom_left     | x,y        | yes      | Bottom-left corner of the channel before extrusion                   |
   * |  1 | top_right       | x,y        | yes      | Top-right corner of the channel before extrusion                     |
   * |  2 | center          | x,y        | yes      | Center of the square obstacle                                        |
   * |  3 | inner_half_side | double     | yes      | Half-side of the inner square obstacle                               |
   * |  4 | outer_half_side | double     | yes      | Half-side of the outer transition square                             |
   * |  5 | rotation_deg    | double     | no       | In-plane rotation (degrees) applied to the inner region (default: 0) |
   * |  6 | pad_bottom      | int        | no       | N cells below the transition (default: 0)                            |
   * |  7 | pad_top         | int        | no       | N cells above the transition (default: 0)                            |
   * |  8 | pad_left        | int        | no       | N cells left of the transition (default: 0)                          |
   * |  9 | pad_right       | int        | no       | N cells right of the transition (default: 0)                         |
   * | 10 | height          | double     | no       | Extrusion height in z, 3D only (default: 1.0)                        |
   * | 11 | n_slices        | int        | no       | Number of z-layers, 3D only (default and min: 2)                     |
   * | 12 | colorize        | true/false | no       | Assign distinct boundary IDs (default: false)                        |
   *
   * Example: @code "0,0 : 10,2 : 5,1 : 0.1 : 0.3 : 15 : 2 : 2 : 5 : 5 : 2. : 2
   * : true "
   * @endcode
   *
   * When @p colorize is true, boundary IDs follow the
   * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   * In 3D, the extruded bottom (z=0) gets id 4 and top (z=height) gets id 5.
   */
  GridUniformChannelWithMeshedSquarePrism(const std::string &grid_arguments);

  /**
   * @brief Generate the 2D/3D channel mesh with a meshed square obstacle.
   *
   * For 2D: Delegates geometry construction to generate_2d_channel_mesh().
   * The boundary IDs are assigned following the subdivided_hyper_rectangle
   * convention when colorize is enabled: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   *
   * For 3D: Delegates geometry construction to generate_2d_channel_mesh(), then
   * extrudes the mesh along the z-axis and attaches the new boundary IDs. The
   * lateral boundary IDs are inherited from the 2D mesh when colorize is
   * enabled. The extruded bottom (z = 0) gets boundary_id = 4 and the top (z =
   * height) gets boundary_id = 5.
   *
   * @param[out] triangulation The triangulation to fill with the channel mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /**
   * @brief Generate the 2D channel mesh geometry with a meshed square obstacle.
   *
   * This private helper creates the composite 2D triangulation by merging:
   * - A modified hyper_ball_balanced from deal.ii so the inner vertex are
   *   placed to form the inner square of half-side @p inner_half_side and the
   *   outer vertices are placed to form the outer square of half-side
   *   @p outer_half_side. When rotated, the inner vertices are rotated by
   *   @p inner_rotation_angle while only half of the outer vertices are rotated
   *   and projected to the closest point on the outer square. The other half
   *   are attached to the outer square corners. The choice of which outer
   *   vertices are rotated and which are attached to the corners is made on a
   *   closest-angle basis to minimize the distortion of the transition cells.
   * - Up to 8 rectangular padding sub-meshes filling the gap between the
   *   transition square and the channel boundary. Those are set up in a way
   *   that their vertices lines up with the vertices projected from the
   *   transition region to avoid the need of introducing new vertices when
   *   merging the triangulations.
   *
   * Material IDs are assigned as follows: 1 for the inner obstacle and 0 for
   * the rest of the channel.
   *
   * If @p colorize is true, boundary IDs are set following the
   * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   *
   * @param[out] triangulation  The triangulation to fill.
   * @param[in]  bottom_left    Bottom-left corner of the channel.
   * @param[in]  top_right      Top-right corner of the channel.
   * @param[in]  center         Center of the square obstacle.
   * @param[in]  inner_half_side   Half-side of the inner square obstacle.
   * @param[in]  outer_half_side   Half-side of the outer transition square.
   * @param[in]  inner_rotation_angle In-plane rotation angle in radians.
   * @param[in]  pad_bottom     Number of subdivisions in the bottom padding.
   * @param[in]  pad_top        Number of subdivisions in the top padding.
   * @param[in]  pad_left       Number of subdivisions in the left padding.
   * @param[in]  pad_right      Number of subdivisions in the right padding.
   * @param[in]  colorize       Whether to assign distinct boundary IDs.
   */
  void
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
                           const bool         colorize);

  /// Point that define the bottom-left of the channel in the xy-plane.
  Point<dim> bottom_left;
  /// Point that define the top-right of the channel in the xy-plane.
  Point<dim> top_right;
  /// Point that define the center of the obstacle in the xy-plane.
  Point<dim> center;
  /// Half-side of the inner square obstacle.
  double inner_half_side;
  /// Half-side of the outer transition square.
  double outer_half_side;
  /// In-plane rotation angle in degrees applied to the inner mesh pattern.
  double inner_rotation_angle;
  /// Number of additional cells to be added to pad the channel from center -
  /// inner_half_side - outer_half_side to the channel bottom boundary (-y
  /// direction).
  unsigned int pad_bottom;
  /// Number of additional cells to be added to pad the channel from center +
  /// inner_half_side + outer_half_side to the channel top boundary (+y
  /// direction).
  unsigned int pad_top;
  /// Number of additional cells to be added to pad the channel from center -
  /// inner_half_side - outer_half_side to the channel left boundary (-x
  /// direction).
  unsigned int pad_left;
  /// Number of additional cells to be added to pad the channel from center +
  /// inner_half_side + outer_half_side to the channel right boundary (+x
  /// direction).
  unsigned int pad_right;
  /// Extrusion height in the z direction, 3D only.
  double height;
  /// Number of layers in the z direction, 3D only. Minimum is 2, which
  /// corresponds to a single layer of cells in the z direction.
  unsigned int n_slices;
  /// Whether to assign distinct boundary IDs to the channel boundaries
  /// following the subdivided_hyper_rectangle convention.
  bool colorize;
};



#endif
