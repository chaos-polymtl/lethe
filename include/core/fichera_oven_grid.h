// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fichera_oven_grid_h
#define lethe_fichera_oven_grid_h

#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <sstream>

using namespace dealii;

/**
 * @brief Generates a 3D Fichera oven mesh by removing cells from a
 * subdivided hyper-rectangle.
 *
 * The geometry is a staircase-shaped domain obtained from a 2x2x3 subdivided
 * hyper-rectangle by removing three groups of cells at increasing z-levels.
 *
 * When @p colorize is true, boundary faces at the top of the remaining
 * staircase (z = top_right[2]) are assigned boundary_id = 1 (e.g., inlet);
 * all other boundary faces retain the default boundary_id = 0.
 */
template <int dim, int spacedim>
class FicheraOvenGrid
{
public:
  /**
   * @brief Constructor that parses geometry parameters from a colon-separated
   * string.
   *
   * @param[in] grid_arguments A colon-separated string, @code
   * "bottom_left:top_right:colorize", with the following fields (coordinates
   * are comma-separated):
   *
   * | #  | Field       | Format     | Required | Description                                   |
   * |----|-------------|------------|----------|-----------------------------------------------|
   * |  0 | bottom_left | x,y,z      | yes      | Bottom-left corner of the bounding box        |
   * |  1 | top_right   | x,y,z      | yes      | Top-right corner of the bounding box          |
   * |  2 | colorize    | true/false | no       | Assign distinct boundary IDs (default: false) |
   *
   * Example: @code "0,0,0 :2,2,3 : true" @endcode
   *
   */
  FicheraOvenGrid(const std::string &grid_arguments);

  /**
   * @brief Generate the Fichera oven mesh by removing cells from a subdivided
   * hyper-rectangle.
   *
   * Starting from a 2x2x3 subdivided hyper-rectangle spanning from
   * @p bottom_left to @p top_right, three groups of cells are removed to create
   * a staircase geometry:
   *
   * @verbatim
   *   Level 1 (z > 1/3 height): remove (+x, +y) quadrant
   *   Level 2 (z > 2/3 height): additionally remove (+x, -y) and (-x, +y)
   * @endverbatim
   *
   * The midpoints along x and y are at (bottom_left + top_right) / 2, and the
   * z-thresholds are at bottom_left[2] + height/3 and bottom_left[2] +
   * 2*height/3.
   *
   * When @p colorize is true, boundary faces at z = top_right[2] are assigned
   * boundary_id = 1 (e.g., inlet/outlet); all other faces keep the default
   * boundary_id = 0.
   *
   * @param[out] triangulation The triangulation to fill with the oven mesh.
   */
  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  Point<dim>  bottom_left;
  Point<dim>  top_right;
  bool        colorize;
};


#endif
