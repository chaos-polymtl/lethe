// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/fichera_oven_grid.h>

template <int dim, int spacedim>
FicheraOvenGrid<dim, spacedim>::FicheraOvenGrid(
  const std::string &grid_arguments)
{
  if constexpr (dim != 3 || spacedim != 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "The Fichera oven mesh is only supported in 3D space with 3D elements."));
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

  // Colorize the inlet. For this geometry, the top face of the remaining
  // staircase (z = top_right[2]) is chosen as the inlet and assigned
  // boundary_id = 1 for the waveguide problem the Fichera oven is used for; all
  // other faces keep the default boundary_id = 0. The top face is identified by
  // checking if the z-coordinate of the face center is close to top_right[2].
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
      "FicheraOvenGrid is only supported for dim = 3 ,spacedim = 3 specializations."));
}

// Explicit template instantiations
template class FicheraOvenGrid<2, 2>;
template class FicheraOvenGrid<2, 3>;
template class FicheraOvenGrid<3, 3>;
