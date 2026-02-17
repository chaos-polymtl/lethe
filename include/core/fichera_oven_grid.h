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
 * @brief Class that creates a custom fichera oven geometry.
 */

template <int dim, int spacedim>
class FicheraOvenGrid
{
public:
  FicheraOvenGrid(const std::string &grid_type, const std::string &grid_arguments);

  void
  make_grid(Triangulation<dim, spacedim> &triangulation);

private:
  std::string grid_arguments;
  Point<dim>   bottom_left;
  Point<dim>   top_right;
  bool         colorize;

};


/**
 * @brief Constructor for the FicheraOvenGrid.
 *
 * @param grid_arguments. A string with 3 parameters
 * @param bottom_left : Point of the lower left corner of the domain
 * @param top_right : Point of the top right corner of the domain
 * @param colorize : Whether to colorize the mesh or not (true/false)
 */

template <int dim, int spacedim>
FicheraOvenGrid<dim, spacedim>::FicheraOvenGrid(const std::string &grid_type,
                                          const std::string &grid_arguments)
{
  if constexpr (dim != 3 || spacedim != 3)
    {
      throw std::runtime_error(
        "The Fichera oven mesh is only supported in 3d space with 3d elements.");
    }

  // Separate arguments of the string (points are given as x,y,z and colorize is given as true/false)
  std::vector<std::string> arguments;
  std::stringstream        s_stream(grid_arguments);
  while (s_stream.good())
    {
      std::string argument;
      getline(s_stream, argument, ':');
      arguments.push_back(argument);
    }

  if (arguments.size() != 3)
    {
      throw std::runtime_error("The Fichera oven grid requires exactly 3 arguments. The required arguments are (lower left point : top right point : colorize). The points should be given as x,y,z and the colorize argument should be given as true/false.");
    }

  // Parse bottom_left point
  std::stringstream bottom_left_stream(arguments[0]);
  std::vector<double> bottom_left_coords;
  std::string       coord_str;
  while (getline(bottom_left_stream, coord_str, ','))
    {
      bottom_left_coords.push_back(std::stod(coord_str));
    }
  if (bottom_left_coords.size() != dim)
    {
      throw std::runtime_error("The bottom left point must have exactly " +
                               std::to_string(dim) + " coordinates using the x,y,z format.");
    }
  for (int i = 0; i < dim; ++i)
    {
      bottom_left[i] = bottom_left_coords[i];
    }

  // Parse top_right point
  std::stringstream top_right_stream(arguments[1]);
  std::vector<double> top_right_coords;
  while (getline(top_right_stream, coord_str, ','))
    {
      top_right_coords.push_back(std::stod(coord_str));
    }
  if (top_right_coords.size() != dim)
    {
      throw std::runtime_error("The top right point must have exactly " +
                               std::to_string(dim) + " coordinates using the x,y,z format.");
    }
  for (int i = 0; i < dim; ++i)
    {
      top_right[i] = top_right_coords[i];
    }

  // Parse colorize flag
  if (arguments[2] == "true")
    {
      colorize = true;
    }
  else if (arguments[2] == "false")
    {
      colorize = false;
    }
  else
    {
      throw std::runtime_error("The colorize argument must be either 'true' or 'false'.");
    }

}

/**
 * @brief make_grid. The make_grid function generates a hyper rectangle of the size of the domain
 * and then transforms it to the Fichera oven geometry. It also colorizes the inlet boundary if specified in the arguments.
 *
 * @param triangulation. The triangulation object on which the grid is generated
 */
template <int dim, int spacedim>
void
FicheraOvenGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> &triangulation)
{
 // Create the domain to be cut
  Triangulation<3> bulk_tria;
  std::vector<unsigned int> repetitions = {2, 2, 3};
  GridGenerator::subdivided_hyper_rectangle(bulk_tria, repetitions, bottom_left,
                                            top_right);
  // Get geometry size in each dimension
  Tensor<1, 3> dimensions = top_right - bottom_left;

  // Remove the undesired cells
  std::set<Triangulation<3>::active_cell_iterator> cells_to_remove;
  for (const auto &cell : bulk_tria.active_cell_iterators()) {
    const Point<3> cell_center = cell->center();

    if ((cell_center[0] > dimensions[0] / 2) &&
        (cell_center[1] > dimensions[1] / 2) &&
        (cell_center[2] > dimensions[2] / 3)) {
      cells_to_remove.insert(cell);
    }
    if ((cell_center[0] > dimensions[0] / 2) &&
        (cell_center[1] < dimensions[1] / 2) &&
        (cell_center[2] > 2. * dimensions[2] / 3)) {
      cells_to_remove.insert(cell);
    }
    if ((cell_center[0] < dimensions[0] / 2) &&
        (cell_center[1] > dimensions[1] / 2) &&
        (cell_center[2] > 2. * dimensions[2] / 3)) {
      cells_to_remove.insert(cell);
    }
  }

  // Create the oven
  GridGenerator::create_triangulation_with_removed_cells(
      bulk_tria, cells_to_remove, triangulation);

  // Colorize the inlet
  if (colorize) {
    for (const auto &face : triangulation.active_face_iterators())
      if (face->at_boundary()) {
        const Point<3> face_center = face->center();
        if (std::abs(face_center[2] - dimensions[2]) < 1e-12)
          face->set_boundary_id(1);
      }
  }
}

#endif
