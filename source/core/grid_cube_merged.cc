// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_cube_merged.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

template <int dim, int spacedim>
GridCubeMerged<dim, spacedim>::GridCubeMerged(const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "Custom cube merged mesh is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;

  // Separate arguments of the string
  std::vector<std::string> arguments;
  std::stringstream        s_stream(grid_arguments);
  while (s_stream.good())
    {
      std::string substr;
      getline(s_stream, substr, ':');
      arguments.push_back(substr);
    }

  // Arguments declaration
  if (arguments.size() != 4)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory cylinder parameters are (plate file name : cylinder file name : cylindrical BID plate : cylindrical BID cylinder)"));
    }
  else
    {
      this->plate_file_name          = arguments[0];
      this->cylinder_file_name       = arguments[1];
      this->cylindrical_bid_plate    = Utilities::string_to_int(arguments[2]);
      this->cylindrical_bid_cylinder = Utilities::string_to_int(arguments[3]);
    }
}


template <>
void
GridCubeMerged<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  // Create a temporary triangulation for the extruded plate with hole
  Triangulation<3, 3> tria_plate;
  GridIn<3, 3>        grid_plate;
  grid_plate.attach_triangulation(tria_plate);
  std::ifstream plate_input_file(this->plate_file_name);
  grid_plate.read_msh(plate_input_file);

  // Create a temporary triangulation for the inner cylinder
  Triangulation<3, 3> tria_cylinder;
  GridIn<3, 3>        grid_cylinder;
  grid_cylinder.attach_triangulation(tria_cylinder);
  std::ifstream cylinder_input_file(this->cylinder_file_name);
  grid_cylinder.read_msh(cylinder_input_file);

  // Set cylindrical manifolds in plate triangulation
  tria_plate.set_manifold(this->cylindrical_bid_plate,
                          CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                 Point<3>({0, 0, 0})));

  for (const auto &cell : tria_plate.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() &&
          face->boundary_id() == this->cylindrical_bid_plate)
        {
          face->set_all_manifold_ids(this->cylindrical_bid_plate);
          cell->set_manifold_id(this->cylindrical_bid_plate);
        }

  // Set cylindrical manifolds in cylinder triangulation
  tria_cylinder.set_manifold(this->cylindrical_bid_cylinder,
                             CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                    Point<3>({0, 0, 0})));

  for (const auto &cell : tria_cylinder.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary() &&
          face->boundary_id() == this->cylindrical_bid_cylinder)
        {
          face->set_all_manifold_ids(this->cylindrical_bid_cylinder);
          cell->set_manifold_id(this->cylindrical_bid_cylinder);
        }

  // Merge triangulations
  GridGenerator::merge_triangulations(
    tria_plate, tria_cylinder, triangulation, 1e-8, true, false);

  // Re-attach cylindrical manifolds
  triangulation.set_manifold(this->cylindrical_bid_plate,
                             CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                    Point<3>({0, 0, 0})));
  triangulation.set_manifold(this->cylindrical_bid_cylinder,
                             CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                    Point<3>({0, 0, 0})));
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
GridCubeMerged<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "GridCubeMerged is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class GridCubeMerged<2, 2>;
template class GridCubeMerged<2, 3>;
template class GridCubeMerged<3, 3>;
