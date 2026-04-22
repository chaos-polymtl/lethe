// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_uniform_channel_with_meshed_cylinder_h
#define lethe_grid_uniform_channel_with_meshed_cylinder_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <sstream>

using namespace dealii;

/**
 * @brief Generates a rectangular channel mesh with a meshed cylinder
 * obstacle. By default, this geometry defines two material IDs. The
 * obstacle is assigned the material_id = 1 and the remaining domain is
 * assigned material_id = 0.
 *
 * The mesh consists of three concentric regions around the cylinder:
 * -# An inner balanced disk of radius @p inner_radius (hyper_ball_balanced).
 * -# A transition shell between @p inner_radius and @p outer_radius, whose
 *    diagonal vertices are adjusted to form a square.
 * -# Rectangular padding regions that fill the gap between the transition
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
 * The manifold objects attached for proper mesh refinement are:
 * - PolarManifold (2D) or CylindricalManifold (3D) on the inner cylinder
 *   boundary faces (manifold_id = 1).
 * - TransfiniteInterpolationManifold on the remaining domain (manifold_id = 0).
 */
template <int dim, int spacedim>
class GridUniformChannelWithMeshedCylinder
{
public:
  /**
   * @brief Constructor that parses geometry parameters from a colon-separated
   * string.
   *
   * @param[in] grid_arguments A colon-separated string with the following
   * fields (coordinates are comma-separated):
   *
   * | #  | Field        | Format     | Required | Description                                        |
   * |----|--------------|------------|----------|----------------------------------------------------|
   * |  0 | bottom_left  | x,y        | yes      | Bottom-left corner of the channel before extrusion |
   * |  1 | top_right    | x,y        | yes      | Top-right corner of the channel before extrusion   |
   * |  2 | center       | x,y        | yes      | Center of the cylinder                             |
   * |  3 | inner_radius | double     | yes      | Radius of the cylinder                             |
   * |  4 | outer_radius | double     | yes      | Radius of the transition region                    |
   * |  5 | pad_bottom   | int        | no       | N cells below the transition (default: 0)          |
   * |  6 | pad_top      | int        | no       | N cells above the transition (default: 0)          |
   * |  7 | pad_left     | int        | no       | N cells left of the transition (default: 0)        |
   * |  8 | pad_right    | int        | no       | N cells right of the transition (default: 0)       |
   * |  9 | height       | double     | no       | Extrusion height in z, 3D only (default: 1.0)      |
   * | 10 | n_slices     | int        | no       | Number of z-layers, 3D only (default and min: 2)   |
   * | 11 | colorize     | true/false | no       | Assign distinct boundary IDs (default: false)      |
   *
   * Example: @code "0,0 : 10,2 : 5,1 : 0.1 : 0.3 : 2 : 2 : 5 : 5 : 2. : 2 :
   * true"
   * @endcode
   *
   * When @p colorize is true, boundary IDs follow the
   * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   * In 3D, the extruded bottom (z=0) gets id 4 and top (z=height) gets id 5.
   */
  GridUniformChannelWithMeshedCylinder(const std::string &grid_arguments);

  /**
   * @brief Generate the 2D or 3D channel mesh with a meshed cylinder and attach
   * manifold objects for proper mesh refinement.
   *
   * For 2D: Delegates geometry construction to generate_2d_channel_mesh(), then
   * attaches:
   * - PolarManifold on the inner-cylinder boundary faces (manifold_id = 1).
   * - TransfiniteInterpolationManifold on the remaining domain (manifold_id =
   * 0).
   * The boundary IDs are assigned following the subdivided_hyper_rectangle
   * convention when colorize is enabled: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   *
   * For 3D: Delegates geometry construction to generate_2d_channel_mesh(), then
   * extrudes the mesh along the z-axis and attaches:
   * - CylindricalManifold on the cylinder surface (manifold_id = 1).
   * - TransfiniteInterpolationManifold on the remaining domain (manifold_id =
   * 0). The lateral boundary IDs are inherited from the 2D mesh when colorize
   * is enabled. The extruded bottom (z = 0) gets boundary_id = 4 and the top (z
   * = height) gets boundary_id = 5.
   *
   * @param[out] triangulation The triangulation to fill with the channel mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /**
   * @brief Generate the 2D channel mesh geometry with a meshed cylinder.
   *
   * This private helper creates the composite 2D triangulation by merging:
   * - A balanced disk of radius @p inner_radius (the cylinder cross-section).
   * - A hyper_shell transition region with squared diagonal vertices between
   *   @p inner_radius and @p outer_radius.
   * - Up to 8 rectangular padding sub-meshes filling the gap between the
   *   transition square and the channel boundary.
   *
   * Manifold integer IDs are assigned (0 = TFI, 1 = polar) but no manifold
   * objects are attached — the caller must do that after this function returns.
   *
   * If @p colorize is true, boundary IDs are set following the
   * subdivided_hyper_rectangle convention: 0 (-x), 1 (+x), 2 (-y), 3 (+y).
   *
   * @param[out] triangulation  The triangulation to fill.
   * @param[in]  bottom_left    Bottom-left corner of the channel.
   * @param[in]  top_right      Top-right corner of the channel.
   * @param[in]  center         Center of the cylinder.
   * @param[in]  inner_radius   Radius of the inner cylinder.
   * @param[in]  outer_radius   Radius of the transition region.
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
                           const double       inner_radius,
                           const double       outer_radius,
                           const unsigned int pad_bottom,
                           const unsigned int pad_top,
                           const unsigned int pad_left,
                           const unsigned int pad_right,
                           const bool         colorize);

  /// Point that define the bottom-left of the channel in the xy-plane.
  Point<dim> bottom_left;
  /// Point that define the top-right of the channel in the xy-plane.
  Point<dim> top_right;
  /// Point that define the center of the cylinder in the xy-plane.
  Point<dim> center;
  /// Radius of the cylinder.
  double inner_radius;
  /// Radius of the transition region between the cylinder and the channel.
  double outer_radius;
  /// Number of additional cells to be added to pad the channel from center -
  /// inner_radius - outer_radius to the channel bottom boundary (-y direction).
  unsigned int pad_bottom;
  /// Number of additional cells to be added to pad the channel from center +
  /// inner_radius + outer_radius to the channel top boundary (+y direction).
  unsigned int pad_top;
  /// Number of additional cells to be added to pad the channel from center -
  /// inner_radius - outer_radius to the channel left boundary (-x direction).
  unsigned int pad_left;
  /// Number of additional cells to be added to pad the channel from center +
  /// inner_radius + outer_radius to the channel right boundary (+x direction).
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
