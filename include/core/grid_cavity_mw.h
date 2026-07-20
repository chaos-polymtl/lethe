// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_grid_cavity_mw_h
#define lethe_grid_cavity_mw_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <array>
#include <sstream>

using namespace dealii;

/** 
* @brief Generates a 3D mesh of a waveguide cavity embedded in a cylinder in which a liquid flows
*
* The geometry consists of three sections joined along the z-axis:
* - A bottom cylinder
* - A central cylinder into which the rectangular-based waveguide is embedded
* - A top cylinder
*
* Boundary IDs are assigned as follows:
* -0: Microwaves (MW) inlet 
* -1: Fluid flow inlet (bottom surface of the cylinder)
* -2: Fluid flow outlet (top surface of the cylinder)
* -3: Walls
*
* Material IDs are assigned as follows:
* -0: Air in the waveguide (paving stone)
* -1: Fluid in the channel (cylinder)
*
 * A CylindricalManifold is attached to the curved lateral surfaces of the cylinder
 * so that mesh refinement preserves the  the shape of the fluid chennel. 
 *
 * @tparam dim The dimension of the mesh, must be 3.
 * @tparam spacedim The dimension of the space, must be 3.
*/

template <int dim, int spacedim>
class GridCavityMw

{
public:
  /**
   * @brief Constructor that parses grid parameters from a colon-separated
   * string.
   *
   * @param[in] grid_arguments A string with the following format:
   * @code "rectangle_width : rectangle_length : outer_radius : inner_length : bottom_height : center_height : top_height : shape" @endcode
   *
   * | #  | Field             | Format        | Required | Description                                                                             |
   * |----|-------------------|---------------|----------|-----------------------------------------------------------------------------------------|
   * |  0 | rectangle_width   | double        | yes      | Waveguide width                                                                         |
   * |  1 | rectangle_length  | double        | yes      | Waveguide length                                                                        |
   * |  2 | outer_radius      | double        | yes      | Radius of the main cylinder                                                             |
   * |  3 | inner_length      | double        | yes      | Depending on the inner shape, it can be the square side or the circle radius            |
   * |  4 | bottom_height     | double        | yes      | Height of the bottom cylinder                                                           |
   * |  5 | center_height     | double        | yes      | Height of the waveguide                                                                 |
   * |  6 | top_height        | double        | yes      | Height of the top cylinder                                                              |
   * |  7 | shape             | true/false    | no       | True (default): The inner shape is a circle False: Inner square                         |
   * Example: @code "2 : 5 : 5 : 3 : 2 : 2 : 3 : true" @endcode
   */
  GridCavityMw(const std::string &grid_arguments);

  /**
   * @brief Generate a cylindrical fluid conduit in which a block-shaped waveguide is embedded.
   *
   * The three parts of the geometry (bottom cylinder, center_cylinder (with the paving stone), top cylinder)
   * are first merged and returned coarse
   *
   * @param[out] triangulation The triangulation to fill with the mesh.
   */

  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  /// Arguments used to generate the grid
  std::string grid_arguments;

  /// Width of the rectangular waveguide
  double rectangle_width;

  /// Length of the rectangular waveguide
  double rectangle_length;

  /// Outer radius of the cylindrical cavity
  double outer_radius;

  /// Inner length (radius if circle, side if square) of the central shape
  double inner_length;

  /// Height of the bottom cylinder
  double bottom_height;

  /// Height of the center section (waveguide embedded in cylinder)
  double center_height;

  /// Height of the top cylinder
  double top_height;

  /// Shape of the inner region: true = circle, false = square
  bool shape;
};


#endif