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
  FicheraOvenGrid(const std::string &grid_arguments);

  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  Point<dim>  bottom_left;
  Point<dim>  top_right;
  bool        colorize;
};


/**
 * @brief Constructor that parses geometry parameters from a colon-separated
 * string.
 *
 * @param[in] grid_arguments A colon-separated string with the following fields
 * (coordinates are comma-separated):
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
template <int dim, int spacedim>
FicheraOvenGrid<dim, spacedim>::FicheraOvenGrid(
  const std::string &grid_arguments)
{
  if constexpr (dim != 3 || spacedim != 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The Fichera oven mesh is only supported in 3d space with 3d elements."));
      return;
    }

  this->grid_arguments = grid_arguments;
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  if (arguments.size() < 2)
    {
      AssertThrow(false,
                  ExcMessage(
                    "The Fichera oven grid requires at least 2 arguments: "
                    "(bottom_left : top_right [: colorize]). "
                    "Points are given as x,y,z and colorize as true/false."));
    }

  // Parse bottom_left point
  std::stringstream   bottom_left_stream(arguments[0]);
  std::vector<double> bottom_left_coords;
  std::string         coord_str;
  while (getline(bottom_left_stream, coord_str, ','))
    bottom_left_coords.push_back(std::stod(coord_str));

  AssertThrow(bottom_left_coords.size() == static_cast<unsigned int>(dim),
              ExcMessage("The bottom_left point must have exactly " +
                         std::to_string(dim) + " coordinates (x,y,z format)."));

  this->bottom_left = Point<dim>(bottom_left_coords[0],
                                 bottom_left_coords[1],
                                 bottom_left_coords[2]);

  // Parse top_right point
  std::stringstream   top_right_stream(arguments[1]);
  std::vector<double> top_right_coords;
  while (getline(top_right_stream, coord_str, ','))
    top_right_coords.push_back(std::stod(coord_str));

  AssertThrow(top_right_coords.size() == static_cast<unsigned int>(dim),
              ExcMessage("The top_right point must have exactly " +
                         std::to_string(dim) + " coordinates (x,y,z format)."));
  this->top_right =
    Point<dim>(top_right_coords[0], top_right_coords[1], top_right_coords[2]);

  // Parse optional colorize flag (defaults to false)
  this->colorize =
    (arguments.size() > 2 && arguments[2] == "true") ? true : false;
}

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
template <>
void
FicheraOvenGrid<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  // Create the initial 2x2x3 subdivided rectangle to be carved
  Triangulation<3>          bulk_tria;
  std::vector<unsigned int> repetitions = {2, 2, 3};
  GridGenerator::subdivided_hyper_rectangle(bulk_tria,
                                            repetitions,
                                            bottom_left,
                                            top_right);

  // Compute the dimensions of the box for later use in cell removal
  const Tensor<1, 3> dimensions = top_right - bottom_left;

  // Remove cells to create the staircase geometry
  std::set<Triangulation<3>::active_cell_iterator> cells_to_remove;
  for (const auto &cell : bulk_tria.active_cell_iterators())
    {
      const Point<3> cell_center = cell->center();

      // Remove cells in the (+x, +y) quadrant above z > 1/3 height
      if ((cell_center[0] > dimensions[0] / 2) &&
          (cell_center[1] > dimensions[1] / 2) &&
          (cell_center[2] > dimensions[2] / 3))
        {
          cells_to_remove.insert(cell);
        }
      // Remove cells in the (+x, -y) above z > 2/3 height
      if ((cell_center[0] > dimensions[0] / 2) &&
          (cell_center[1] < dimensions[1] / 2) &&
          (cell_center[2] > 2. * dimensions[2] / 3))
        {
          cells_to_remove.insert(cell);
        }
      // Remove cells in the (-x, +y) above z > 2/3 height
      if ((cell_center[0] < dimensions[0] / 2) &&
          (cell_center[1] > dimensions[1] / 2) &&
          (cell_center[2] > 2. * dimensions[2] / 3))
        {
          cells_to_remove.insert(cell);
        }
    }

  // Create the oven
  GridGenerator::create_triangulation_with_removed_cells(bulk_tria,
                                                         cells_to_remove,
                                                         triangulation);

  // Colorize the inlet
  if (colorize)
    {
      for (const auto &face : triangulation.active_face_iterators())
        if (face->at_boundary())
          {
            const Point<3> face_center = face->center();
            if (std::abs(face_center[2] - dimensions[2]) < 1e-12)
              face->set_boundary_id(1);
          }
    }
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
FicheraOvenGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "FicheraOvenGrid is only supported for <3,3> <dim,spacedim> specializations."));
}

#endif
