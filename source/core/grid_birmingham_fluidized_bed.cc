// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grid_birmingham_fluidized_bed.h>

template <int dim, int spacedim>
BirminghamFluidizedBedGrid<dim, spacedim>::BirminghamFluidizedBedGrid(
  const std::string &grid_arguments)
{
  if constexpr (!(dim == 3 && spacedim == 3))
    {
      AssertThrow(
        false,
        ExcMessage(
          "The Birmingham fluidized bed mesh is only supported in 3D space with 3D elements."));
      return;
    }

  this->grid_arguments = grid_arguments;
}


template <>
void
BirminghamFluidizedBedGrid<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  const double       r_small       = 0.154 / 2;
  const double       r_large       = 0.254 / 2;
  const double       half_len_bot  = 0.924 / 2;
  const double       half_len_cone = 0.285 / 2;
  const double       half_len_top  = 0.5 / 2;
  const unsigned int n_axial_bot   = 10;
  const unsigned int n_axial_top   = 5;

  // Bottom cylinder (small radius), shifted to [0, 2*half_len_bot]
  Triangulation<3> tria_bottom;
  GridGenerator::subdivided_cylinder(tria_bottom,
                                     n_axial_bot,
                                     r_small,
                                     half_len_bot);
  GridTools::shift(Point<3>(half_len_bot, 0, 0), tria_bottom);

  // Truncated cone from r_small to r_large,
  // shifted to [2*half_len_bot, 2*half_len_bot + 2*half_len_cone]
  Triangulation<3> tria_cone;
  GridGenerator::truncated_cone(tria_cone, r_small, r_large, half_len_cone);
  GridTools::shift(Point<3>(2.0 * half_len_bot + half_len_cone, 0, 0),
                   tria_cone);

  // Top cylinder (large radius = 2x small),
  // shifted to [2*half_len_bot + 2*half_len_cone, ...]
  Triangulation<3> tria_top;
  GridGenerator::subdivided_cylinder(tria_top,
                                     n_axial_top,
                                     r_large,
                                     half_len_top);
  GridTools::shift(
    Point<3>(2.0 * half_len_bot + 2.0 * half_len_cone + half_len_top, 0, 0),
    tria_top);

  // Merge: bottom + cone -> temp, then temp + top -> result
  Triangulation<3> tria_temp;
  GridGenerator::merge_triangulations(tria_bottom, tria_cone, tria_temp);
  GridGenerator::merge_triangulations(tria_temp, tria_top, triangulation);

  // Re-assign manifold IDs lost during merge.
  // End-cap faces (all vertices at x ≈ 0 or x ≈ total_length) stay flat;
  // all other boundary faces/edges are curved and get manifold_id 1.
  const double total_length =
    2.0 * (half_len_bot + half_len_cone + half_len_top);
  const double tol = 1e-10;

  for (auto &cell : triangulation.active_cell_iterators())
    {
      cell->set_manifold_id(0);
      for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
        {
          if (!cell->face(f)->at_boundary())
            continue;

          bool is_end_cap = true;
          for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face; ++v)
            {
              const double x = cell->face(f)->vertex(v)[0];
              if (std::abs(x) > tol && std::abs(x - total_length) > tol)
                {
                  is_end_cap = false;
                  break;
                }
            }

          if (!is_end_cap)
            {
              cell->face(f)->set_manifold_id(1);
              for (unsigned int l = 0; l < GeometryInfo<3>::lines_per_face; ++l)
                cell->face(f)->line(l)->set_manifold_id(1);
              cell->face(f)->set_boundary_id(0);
            }
          else
            {
              const double x = cell->face(f)->center()[0];
              if (std::abs(x) < tol)
                cell->face(f)->set_boundary_id(1);
              else
                cell->face(f)->set_boundary_id(2);
            }
        }
    }

  // CylindricalManifold along the x-axis handles both the cylinder
  // surfaces and the cone (whose radius varies linearly, matching the
  // manifold's linear radius averaging).
  triangulation.set_manifold(1, CylindricalManifold<3>(0));
  triangulation.set_manifold(0, FlatManifold<3>());
}

// Fallback make_grid definition for unsupported template parameters. This
// provides a linker-visible symbol and a clear runtime error when the
// class is instantiated for dim/spacedim combinations that are not
// specialized above.
template <int dim, int spacedim>
void
BirminghamFluidizedBedGrid<dim, spacedim>::make_grid(
  Triangulation<dim, spacedim> & /*triangulation*/)
{
  AssertThrow(
    false,
    ExcMessage(
      "BirminghamFluidizedBedGrid is only implemented for dim = 3 and spacedim = 3."));
}

// Explicit template instantiations
template class BirminghamFluidizedBedGrid<2, 2>;
template class BirminghamFluidizedBedGrid<2, 3>;
template class BirminghamFluidizedBedGrid<3, 3>;
