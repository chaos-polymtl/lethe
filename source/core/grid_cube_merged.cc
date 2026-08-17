// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_cube_merged.h>
#include <core/lethe_grid_tools.h>

#include <deal.II/grid/grid_generator.h>
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
  std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  // Arguments declaration
  if (arguments.size() != 4)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory cylinder parameters are (R_cylinder : L_plate : H : n_sudivisions)"));
    }
  else
    {
      this->cylinder_radius   = Utilities::string_to_double(arguments[0]);
      this->plate_half_length = Utilities::string_to_double(arguments[1]);
      this->height            = Utilities::string_to_double(arguments[2]);
      this->n_subdivisions    = Utilities::string_to_int(arguments[3]);

      // To match the number of cells at the cylindrical interface, the cylinder
      // geometry is refined once before being merged into the extruded plate.
      // The prescribed number of subdivisions, therefore, needs to be divided
      // by half to create the initial cylinder. Hence, he number of
      // subdivisions needs to be an even number
      AssertThrow(this->n_subdivisions % 2 == 0,
                  ExcMessage(
                    "The number of subdivisions needs to be an even number."));
    }
}


template <>
void
GridCubeMerged<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  // Create grid of extruded plate with hole
  Triangulation<3, 3> tria_plate;
  GridGenerator::hyper_cube_with_cylindrical_hole(tria_plate,
                                                  cylinder_radius,
                                                  plate_half_length,
                                                  height,
                                                  n_subdivisions,
                                                  false);

  // Create grid of cylinder
  Triangulation<3, 3> tria_cylinder;
  GridGenerator::subdivided_cylinder(tria_cylinder,
                                     n_subdivisions / 2,
                                     cylinder_radius,
                                     height / 2);

  // Rotate the cylinder so that it is aligned with the z axis
  GridTools::rotate(Tensor<1, 3>({0, 1, 0}),
                    std::numbers::pi / 2,
                    tria_cylinder);
  // Shift the cylinder so that its 2D extrusion plane is the same as tria_plate
  GridTools::shift(Tensor<1, 3>({0, 0, 1}), tria_cylinder);
  // Remove cylindrical manifold in x direction
  tria_cylinder.reset_all_manifolds();
  // Set new manifold on z direction (manifold ID for the hull of the cylinder
  // is 0)
  tria_cylinder.set_manifold(0,
                             CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                    Point<3>()));

  // Refine the cylinder mesh once so that the number of cells at the interface
  // match with tria_plate
  tria_cylinder.refine_global(1);

  // Flatten cylinder triangulation so that it can be merged with the extruded
  // plate with hole
  Triangulation<3, 3> tria_cylinder_flat;
  GridGenerator::flatten_triangulation(tria_cylinder, tria_cylinder_flat);

  // Merge triangulations
  GridGenerator::merge_triangulations(tria_plate,
                                      tria_cylinder_flat,
                                      triangulation,
                                      1e-8);

  // Re-assign cylindrical manifold
  double tol = 1e-10;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
        {
          bool is_cylindrical_boundary = true;
          for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; ++v)
            {
              const auto vertex = cell->face(f)->vertex(v);

              // Check if vertex lies on the former cylindrical boundary
              double distance = LetheGridTools::find_point_line_distance(
                Point<3>(), Tensor<1, 3>({0, 0, 1}), vertex);
              if (std::abs(distance - this->cylinder_radius) > tol)
                {
                  is_cylindrical_boundary = false;
                  break;
                }
            }

          if (is_cylindrical_boundary)
            {
              cell->face(f)->set_manifold_id(1);
              for (unsigned int l = 0; l < GeometryInfo<3>::lines_per_face; ++l)
                cell->face(f)->line(l)->set_manifold_id(1);
            }
        }
    }

  triangulation.set_manifold(0, FlatManifold<3>());
  triangulation.set_manifold(1,
                             CylindricalManifold<3>(Tensor<1, 3>({0, 0, 1}),
                                                    Point<3>()));
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
