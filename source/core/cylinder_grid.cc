// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/cylinder_grid.h>

template <int dim, int spacedim>
CylinderGrid<dim, spacedim>::CylinderGrid(const std::string &grid_type,
                                          const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "Custom cylinder mesh is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;
  if (grid_type == "cylinder_classic")
    this->cylinder_type = CylinderType::classic;
  else if (grid_type == "cylinder_balanced")
    this->cylinder_type = CylinderType::balanced;
  else if (grid_type == "cylinder_squared")
    this->cylinder_type = CylinderType::squared;
  else if (grid_type == "cylinder_regularized")
    this->cylinder_type = CylinderType::regularized;
  else
    AssertThrow(
      false,
      ExcMessage(
        "Unsupported custom cylinder, the choices are: cylinder_classic|cylinder_balanced|cylinder_squared|cylinder_regularized"));


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
  if (arguments.size() != 3)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Mandatory cylinder parameters are (x subdivisions: radius : half height)"));
    }
  else
    {
      std::vector<double> arguments_double =
        Utilities::string_to_double(arguments);
      this->subdivisions = static_cast<int>(arguments_double[0]);
      this->radius       = arguments_double[1];
      this->half_height  = arguments_double[2];
    }
}


template <>
void
CylinderGrid<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  if (cylinder_type == CylinderType::classic)
    {
      // Create a subdivided cylinder from deal.ii
      GridGenerator::subdivided_cylinder(triangulation,
                                         subdivisions,
                                         radius,
                                         half_height);
    }
  else
    {
      // Create a temporary 2D mesh
      Triangulation<2, 2> temporary_triangulation;

      // Create a spherical manifold for 2D mesh
      Point<2>                      center(0.0, 0.0);
      const SphericalManifold<2, 2> m0(center);

      if (cylinder_type == CylinderType::regularized ||
          cylinder_type == CylinderType::squared)
        {
          // Create a square mesh
          double real_radius = radius * std::sin(M_PI_4);
          GridGenerator::hyper_cube(temporary_triangulation,
                                    -real_radius,
                                    real_radius,
                                    true);

          // Assign boundary 0 to perimeter as for cylinder
          for (const auto &cell :
               temporary_triangulation.active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  // Looping through all the faces of the cell
                  for (const auto &face : cell->face_iterators())
                    {
                      // Check to see if the face is located at boundary
                      if (face->at_boundary())
                        {
                          face->set_boundary_id(0);
                        }
                    }
                }
            }
        }
      else if (cylinder_type == CylinderType::balanced)
        {
          GridGenerator::hyper_ball_balanced(temporary_triangulation,
                                             center,
                                             radius);
        }

      temporary_triangulation.reset_all_manifolds();
      temporary_triangulation.set_all_manifold_ids_on_boundary(0);
      temporary_triangulation.set_manifold(0, m0);

      if (cylinder_type == CylinderType::regularized)
        {
          // Pre-refinement to reduce mesh size at corners before
          // reinitialization
          temporary_triangulation.refine_global(2);
          GridTools::regularize_corner_cells(temporary_triangulation);

          // Flatten the triangulation
          Triangulation<2, 2> flat_temporary_triangulation;
          flat_temporary_triangulation.copy_triangulation(
            temporary_triangulation);
          temporary_triangulation.clear();
          GridGenerator::flatten_triangulation(flat_temporary_triangulation,
                                               temporary_triangulation);
        }

      // Extrude the 2D temporary mesh to 3D cylinder
      GridGenerator::extrude_triangulation(temporary_triangulation,
                                           subdivisions + 1,
                                           2.0 * half_height,
                                           triangulation,
                                           true);

      // Rotate mesh in x-axis and set the (0,0,0) at the barycenter
      // to be comparable to dealii cylinder meshes
      Tensor<1, 3> axis_vector({0.0, 1.0, 0.0});
      GridTools::rotate(axis_vector, M_PI_2, triangulation);
      Tensor<1, 3> shift_vector({-half_height, 0.0, 0.0});
      GridTools::shift(shift_vector, triangulation);

      // Force the manifold id to be zero in the case of the balanced
      // cylinder
      if (cylinder_type == CylinderType::balanced)
        triangulation.reset_manifold(1);

      // Add a cylindrical manifold on the final unrefined mesh
      const CylindricalManifold<3, 3> m1(0);
      triangulation.set_manifold(0, m1);
    }
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
CylinderGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "CylinderGrid is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class CylinderGrid<2, 2>;
template class CylinderGrid<2, 3>;
template class CylinderGrid<3, 3>;
