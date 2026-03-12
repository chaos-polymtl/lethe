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

  // Parse optional enable_chimney flag (defaults to true)
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');

  this->enable_chimney = !(arguments.size() > 0 && arguments[0] == "false");
}


template <>
void
BirminghamFluidizedBedGrid<3, 3>::make_grid(Triangulation<3, 3> &triangulation)
{
  const double       r_small               = 0.154 / 2;
  const double       r_large               = 0.254 / 2;
  const double       half_len_bot          = 0.924 / 2;
  const double       half_len_cone         = 0.285 / 2;
  const double       half_len_top          = 0.5 / 2;
  const double       chimney_length        = 0.1;
  const unsigned int n_chimney_refinements = 1;
  const unsigned int n_axial_bot           = 10;
  const unsigned int n_axial_top           = 5;
  const unsigned int n_axial_chimney       = 2;
  const double       tol                   = 1e-10;

  // ---- 1. Build the cylindrical geometry ----

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

  // Top cylinder (large radius),
  // shifted to [2*half_len_bot + 2*half_len_cone, ...]
  Triangulation<3> tria_top;
  GridGenerator::subdivided_cylinder(tria_top,
                                     n_axial_top,
                                     r_large,
                                     half_len_top);
  GridTools::shift(
    Point<3>(2.0 * half_len_bot + 2.0 * half_len_cone + half_len_top, 0, 0),
    tria_top);

  // Length of the cylindrical section (before chimney)
  const double base_length =
    2.0 * (half_len_bot + half_len_cone + half_len_top);

  // Merge bottom + cone + top
  Triangulation<3> tria_temp1;
  GridGenerator::merge_triangulations(tria_bottom, tria_cone, tria_temp1);
  Triangulation<3> tria_merged;
  GridGenerator::merge_triangulations(tria_temp1, tria_top, tria_merged);

  if (enable_chimney)
    {
      // ---- 2. Attach manifolds and refine ----

      // Assign manifold IDs to curved surfaces before refining so that new
      // vertices are placed on the cylindrical/conical surface.
      for (auto &cell : tria_merged.active_cell_iterators())
        {
          cell->set_manifold_id(0);
          for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
            {
              if (!cell->face(f)->at_boundary())
                continue;

              bool is_end_cap = true;
              for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face;
                   ++v)
                {
                  const double x = cell->face(f)->vertex(v)[0];
                  if (std::abs(x) > tol && std::abs(x - base_length) > tol)
                    {
                      is_end_cap = false;
                      break;
                    }
                }
              if (!is_end_cap)
                {
                  cell->face(f)->set_manifold_id(1);
                  for (unsigned int l = 0; l < GeometryInfo<3>::lines_per_face;
                       ++l)
                    cell->face(f)->line(l)->set_manifold_id(1);
                }
            }
        }

      tria_merged.set_manifold(1, CylindricalManifold<3>(0));
      tria_merged.set_manifold(0, FlatManifold<3>());

      tria_merged.refine_global(n_chimney_refinements);

      // ---- 3. Flatten the refined mesh ----

      // flatten_triangulation turns every active cell of the refined mesh
      // into a coarse cell, so the result can be merged with the chimney.
      Triangulation<3> tria_flat;
      GridGenerator::flatten_triangulation(tria_merged, tria_flat);

      // ---- 4. Find an off-center face on the top end cap ----

      // Target point in the first shell ring, offset in the +z direction.
      const Point<3> chimney_target(base_length, 0.0, r_large / 2.0);

      double                  best_dist_sq = r_large * r_large + 1.0;
      std::array<Point<3>, 4> best_face_vertices;

      for (const auto &cell : tria_flat.active_cell_iterators())
        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
          {
            if (!cell->face(f)->at_boundary())
              continue;
            const auto &face = cell->face(f);

            // Only consider faces on the top end cap (x ≈ base_length)
            bool on_top_cap = true;
            for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face;
                 ++v)
              if (std::abs(face->vertex(v)[0] - base_length) > tol)
                {
                  on_top_cap = false;
                  break;
                }
            if (!on_top_cap)
              continue;

            const Point<3> center  = face->center();
            const double   dist_sq = (center - chimney_target).norm_square();

            if (dist_sq < best_dist_sq)
              {
                best_dist_sq = dist_sq;
                for (unsigned int v = 0; v < 4; ++v)
                  best_face_vertices[v] = face->vertex(v);
              }
          }

      AssertThrow(best_dist_sq < r_large * r_large,
                  ExcMessage(
                    "Could not find a suitable face on the top end cap "
                    "for the chimney placement."));

      // ---- 5. Build the chimney triangulation ----

      // Sort the four face vertices into (y, z) quadrants relative to the
      // face center so that the hex vertex ordering
      // (index = x_bit + 2*y_bit + 4*z_bit) is satisfied.
      Point<3> face_center;
      for (const auto &p : best_face_vertices)
        face_center += p;
      face_center /= 4.0;

      std::array<Point<3>, 4> sorted_fv;
      for (const auto &p : best_face_vertices)
        {
          unsigned int idx = 0;
          if (p[1] > face_center[1])
            idx += 1; // high y -> y_bit = 1
          if (p[2] > face_center[2])
            idx += 2; // high z -> z_bit = 1
          sorted_fv[idx] = p;
        }

      // Create chimney vertices: (n_axial_chimney + 1) layers of 4
      // vertices, each layer extruded along the x-axis.
      const unsigned int    n_layers = n_axial_chimney;
      std::vector<Point<3>> chimney_vertices(4 * (n_layers + 1));
      for (unsigned int layer = 0; layer <= n_layers; ++layer)
        {
          const double dx = chimney_length * layer / n_layers;
          for (unsigned int v = 0; v < 4; ++v)
            {
              chimney_vertices[layer * 4 + v] = sorted_fv[v];
              chimney_vertices[layer * 4 + v][0] += dx;
            }
        }

      // Create chimney hex cells.  Each cell spans one layer; the vertex
      // layout follows deal.II convention
      // (index = x_bit + 2*y_bit + 4*z_bit).
      std::vector<CellData<3>> chimney_cells(n_layers);
      for (unsigned int layer = 0; layer < n_layers; ++layer)
        {
          const unsigned int L = layer * 4;       // left  layer base index
          const unsigned int R = (layer + 1) * 4; // right layer base index
          chimney_cells[layer].vertices = {L + 0,
                                           R + 0, // hex 0, 1
                                           L + 1,
                                           R + 1, // hex 2, 3
                                           L + 2,
                                           R + 2, // hex 4, 5
                                           L + 3,
                                           R + 3}; // hex 6, 7
        }

      Triangulation<3> tria_chimney;
      tria_chimney.create_triangulation(chimney_vertices,
                                        chimney_cells,
                                        SubCellData());

      // ---- 6. Final merge ----

      GridGenerator::merge_triangulations(tria_flat,
                                          tria_chimney,
                                          triangulation);

      // ---- 7. Reassign manifold and boundary IDs ----

      // Faces with all vertices at x >= base_length are chimney walls or
      // the annular top wall; they stay flat. The inlet end cap (x ≈ 0) is
      // also flat. All other non-end-cap faces are curved and receive
      // manifold_id 1 (CylindricalManifold).
      const double total_length = base_length + chimney_length;

      for (auto &cell : triangulation.active_cell_iterators())
        {
          cell->set_manifold_id(0);
          for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
            {
              if (!cell->face(f)->at_boundary())
                continue;

              bool is_inlet          = true;
              bool is_outlet         = true;
              bool at_or_beyond_base = true;

              for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face;
                   ++v)
                {
                  const double x = cell->face(f)->vertex(v)[0];
                  if (std::abs(x) > tol)
                    is_inlet = false;
                  if (std::abs(x - total_length) > tol)
                    is_outlet = false;
                  if (x < base_length - tol)
                    at_or_beyond_base = false;
                }

              if (is_inlet)
                {
                  cell->face(f)->set_boundary_id(1);
                }
              else if (is_outlet)
                {
                  cell->face(f)->set_boundary_id(2);
                }
              else if (at_or_beyond_base)
                {
                  // Chimney lateral walls or annular top wall around the
                  // chimney opening. These are flat surfaces.
                  cell->face(f)->set_boundary_id(0);
                }
              else
                {
                  // Curved lateral surfaces of cylinders and cone
                  cell->face(f)->set_manifold_id(1);
                  for (unsigned int l = 0; l < GeometryInfo<3>::lines_per_face;
                       ++l)
                    cell->face(f)->line(l)->set_manifold_id(1);
                  cell->face(f)->set_boundary_id(0);
                }
            }
        }

      triangulation.set_manifold(1, CylindricalManifold<3>(0));
      triangulation.set_manifold(0, FlatManifold<3>());
    }
  else
    {
      // ---- No chimney: simple cylindrical geometry ----

      triangulation.copy_triangulation(tria_merged);

      for (auto &cell : triangulation.active_cell_iterators())
        {
          cell->set_manifold_id(0);
          for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f)
            {
              if (!cell->face(f)->at_boundary())
                continue;

              bool is_end_cap = true;
              for (unsigned int v = 0; v < GeometryInfo<3>::vertices_per_face;
                   ++v)
                {
                  const double x = cell->face(f)->vertex(v)[0];
                  if (std::abs(x) > tol && std::abs(x - base_length) > tol)
                    {
                      is_end_cap = false;
                      break;
                    }
                }

              if (!is_end_cap)
                {
                  cell->face(f)->set_manifold_id(1);
                  for (unsigned int l = 0; l < GeometryInfo<3>::lines_per_face;
                       ++l)
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

      triangulation.set_manifold(1, CylindricalManifold<3>(0));
      triangulation.set_manifold(0, FlatManifold<3>());
    }
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
