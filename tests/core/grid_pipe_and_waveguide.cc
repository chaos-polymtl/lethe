// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the generation of the mesh of a waveguide embedded in a
 * cylinder. The grid is only implemented in 3D, so only dim = spacedim = 3 is
 * tested.
 *
 * Two layouts of the cross section are tested: the simplest one (an inner
 * square divided into 2 columns and 3 rows, a single layer of ring cells), and
 * the layout of the heated rotating drum example, with more cells, pulled
 * corners and a band of 3 rows.
 *
 * The number of boundary faces per boundary id is reported on the coarse mesh.
 * The mesh is then refined, and the test reports:
 * - whether all the cells have a positive measure;
 * - the volume of the mesh, next to the exact volume of the union of the
 *   cylinder and of the waveguide, which it approaches as the mesh is refined;
 * - whether the lateral wall of the cylinder, including the faces shared by
 *   the cylinder and the waveguide, lies on the cylinder;
 * - whether the walls of the waveguide are planar;
 * - the number of vertices on the sides of the inner rectangle, which shows
 *   that its refined edges stay straight.
 */

// Deal.II
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_pipe_and_waveguide.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <utility>
#include <vector>

using Grid = GridPipeAndWaveguide<3, 3>;

/**
 * @brief Check that the lateral wall of the cylinder lies on the cylinder.
 *
 * The faces checked are the boundary faces of the lateral wall and the
 * interior faces shared by the cylinder and the waveguide. All their vertices
 * must be at a distance outer_radius from the axis of the cylinder.
 *
 * @param[in] triangulation The mesh to check.
 * @param[in] center Position of the axis of the cylinder in the (x, y) plane.
 * @param[in] outer_radius Radius of the cylinder.
 * @param[in] tolerance Tolerance on the distance to the axis.
 *
 * @return True if all the vertices of these faces are on the cylinder.
 */
bool
lateral_wall_is_on_cylinder(const Triangulation<3, 3> &triangulation,
                            const Point<2>            &center,
                            const double               outer_radius,
                            const double               tolerance)
{
  bool all_vertices_on_cylinder = true;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->material_id() != Grid::cylinder_material_id)
        continue;

      for (const auto f : cell->face_indices())
        {
          const bool face_is_on_lateral_wall =
            cell->face(f)->at_boundary() ?
              cell->face(f)->boundary_id() == Grid::cylinder_wall_boundary_id :
              cell->neighbor(f)->material_id() == Grid::waveguide_material_id;
          if (!face_is_on_lateral_wall)
            continue;

          for (const auto v : cell->face(f)->vertex_indices())
            {
              const Point<3> &vertex = cell->face(f)->vertex(v);
              const double    distance_to_axis =
                center.distance(Point<2>(vertex[0], vertex[1]));
              if (std::abs(distance_to_axis - outer_radius) > tolerance)
                all_vertices_on_cylinder = false;
            }
        }
    }
  return all_vertices_on_cylinder;
}

/**
 * @brief Check that the walls of the waveguide are planar.
 *
 * Each boundary face of the walls of the waveguide must lie on one of the
 * planes y = 0, y = rectangle_width, z = 0 and z = middle_height.
 *
 * @param[in] triangulation The mesh to check.
 * @param[in] rectangle_width Width of the waveguide.
 * @param[in] middle_height Height of the waveguide.
 * @param[in] tolerance Tolerance on the distance to the planes.
 *
 * @return True if each face of the walls of the waveguide lies on one of these
 * planes.
 */
bool
waveguide_walls_are_planar(const Triangulation<3, 3> &triangulation,
                           const double               rectangle_width,
                           const double               middle_height,
                           const double               tolerance)
{
  // Planes of the walls of the waveguide, as (component, position) pairs
  const std::vector<std::pair<unsigned int, double>> wall_planes = {
    {1, 0.}, {1, rectangle_width}, {2, 0.}, {2, middle_height}};

  bool all_faces_planar = true;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto f : cell->face_indices())
      {
        if (!cell->face(f)->at_boundary() ||
            cell->face(f)->boundary_id() != Grid::waveguide_wall_boundary_id)
          continue;

        bool face_is_on_a_wall_plane = false;
        for (const auto &[component, position] : wall_planes)
          {
            bool all_vertices_on_plane = true;
            for (const auto v : cell->face(f)->vertex_indices())
              if (std::abs(cell->face(f)->vertex(v)[component] - position) >
                  tolerance)
                all_vertices_on_plane = false;
            if (all_vertices_on_plane)
              face_is_on_a_wall_plane = true;
          }
        if (!face_is_on_a_wall_plane)
          all_faces_planar = false;
      }
  return all_faces_planar;
}

/**
 * @brief Count the vertices that lie on the sides of the inner rectangle.
 *
 * The inner rectangle [-inner_half_width, inner_half_width] x
 * [-inner_half_height, inner_half_height] is centered on the axis of the
 * cylinder. Its sides stay straight after refinement only if the vertices
 * added on them stay on the lines x = +/- inner_half_width and
 * y = +/- inner_half_height, so that the number of vertices on these lines
 * would drop if they were bent. The corners pulled inward are not on them.
 *
 * @param[in] triangulation The mesh to check.
 * @param[in] center Position of the axis of the cylinder in the (x, y) plane.
 * @param[in] inner_half_width Half-width of the inner rectangle.
 * @param[in] inner_half_height Half-height of the inner rectangle.
 * @param[in] tolerance Tolerance on the distance to the sides.
 *
 * @return The number of vertices on the sides of the inner rectangle.
 */
unsigned int
count_vertices_on_inner_rectangle(const Triangulation<3, 3> &triangulation,
                                  const Point<2>            &center,
                                  const double               inner_half_width,
                                  const double               inner_half_height,
                                  const double               tolerance)
{
  unsigned int n_vertices_on_inner_rectangle = 0;
  for (const auto &vertex : triangulation.get_vertices())
    {
      const double x_from_axis = std::abs(vertex[0] - center[0]);
      const double y_from_axis = std::abs(vertex[1] - center[1]);

      const bool on_left_or_right_side =
        std::abs(x_from_axis - inner_half_width) < tolerance &&
        y_from_axis < inner_half_height + tolerance;
      const bool on_top_or_bottom_side =
        std::abs(y_from_axis - inner_half_height) < tolerance &&
        x_from_axis < inner_half_width + tolerance;
      if (on_left_or_right_side || on_top_or_bottom_side)
        ++n_vertices_on_inner_rectangle;
    }
  return n_vertices_on_inner_rectangle;
}

void
run_test(const std::string &case_name,
         const std::string &grid_arguments,
         const unsigned int n_refinements)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  // Geometry, in the order of the grid arguments
  const std::vector<std::string> arguments =
    Utilities::split_string_list(grid_arguments, ':');
  const double rectangle_width   = Utilities::string_to_double(arguments[0]);
  const double rectangle_length  = Utilities::string_to_double(arguments[1]);
  const double outer_radius      = Utilities::string_to_double(arguments[2]);
  const double inner_half_width  = Utilities::string_to_double(arguments[3]);
  const double inner_half_height = Utilities::string_to_double(arguments[4]);
  const double bottom_height     = Utilities::string_to_double(arguments[11]);
  const double middle_height     = Utilities::string_to_double(arguments[12]);
  const double top_height        = Utilities::string_to_double(arguments[13]);

  Triangulation<3, 3> triangulation;
  Grid                grid(grid_arguments);
  grid.make_grid(triangulation);

  deallog << "Number of coarse cells : " << triangulation.n_active_cells()
          << std::endl;

  // Count the number of faces per boundary id on the coarse mesh
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const auto f : cell->face_indices())
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  triangulation.refine_global(n_refinements);
  deallog << "Number of active cells after " << n_refinements
          << " refinements : " << triangulation.n_active_cells() << std::endl;

  double minimum_cell_measure = std::numeric_limits<double>::max();
  for (const auto &cell : triangulation.active_cell_iterators())
    minimum_cell_measure = std::min(minimum_cell_measure, cell->measure());
  deallog << "All cell measures are positive : "
          << (minimum_cell_measure > 0. ? "yes" : "no") << std::endl;

  /*
  The exact volume is the volume of the cylinder plus the volume of the part of
  the waveguide outside the cylinder. The waveguide ends on the chord
  x = rectangle_length of the circle, at a distance
  distance_center_to_chord = sqrt(outer_radius^2 - rectangle_width^2 / 4) from
  the center. The circular segment between this chord and the circle, of area
  outer_radius^2 acos(distance_center_to_chord / outer_radius)
  - distance_center_to_chord rectangle_width / 2,
  is inside the cylinder and must not be counted twice.
  */
  const double distance_center_to_chord = std::sqrt(
    outer_radius * outer_radius - rectangle_width * rectangle_width / 4.);
  const double segment_area =
    outer_radius * outer_radius *
      std::acos(distance_center_to_chord / outer_radius) -
    distance_center_to_chord * rectangle_width / 2.;
  const double exact_volume =
    std::numbers::pi * outer_radius * outer_radius *
      (bottom_height + middle_height + top_height) +
    (rectangle_length * rectangle_width - segment_area) * middle_height;

  deallog << "Mesh volume  : " << GridTools::volume(triangulation) << std::endl;
  deallog << "Exact volume : " << exact_volume << std::endl;

  // Axis of the cylinder, at the same position as in make_grid
  const Point<2> center(rectangle_length + distance_center_to_chord,
                        rectangle_width / 2.);
  const double   tolerance = 1e-10 * outer_radius;

  deallog << "Lateral wall on the cylinder : "
          << (lateral_wall_is_on_cylinder(
                triangulation, center, outer_radius, tolerance) ?
                "yes" :
                "no")
          << std::endl;
  deallog << "Waveguide walls planar : "
          << (waveguide_walls_are_planar(
                triangulation, rectangle_width, middle_height, tolerance) ?
                "yes" :
                "no")
          << std::endl;
  deallog << "Number of vertices on the sides of the inner rectangle : "
          << count_vertices_on_inner_rectangle(triangulation,
                                               center,
                                               inner_half_width,
                                               inner_half_height,
                                               tolerance)
          << std::endl;

  // Write the refined mesh in VTK format, so that it can be opened in ParaView
  // as a visual check
  GridOut           go;
  const std::string vtk_filename =
    "grid_pipe_and_waveguide_" + case_name + ".vtk";
  std::ofstream vtk_out(vtk_filename);
  go.write_vtk(triangulation, vtk_out);
}

int
main()
{
  try
    {
      initlog();

      // Input arguments:
      // rectangle width   | rectangle length  | outer radius
      // inner half-width  | inner half-height | corner angle
      // corner pull       | n columns         | n strip rows
      // n band rows       | n ring layers     | bottom height
      // middle height     | top height        | bottom resolution
      // middle resolution | top resolution
      run_test("simple_layout",
               "0.07 : 0.05 : 0.1 : 0.05 : 0.05 : 45 : 0 : 2 : 1 : 1 : 1 : "
               "0.1 : 0.07 : 0.1 : 3 : 2 : 3",
               2);
      run_test("heated_drum_layout",
               "0.04318 : 0.1 : 0.069 : 0.0363 : 0.0363 : 39.2 : 0.27 : 5 : "
               "1 : 3 : 2 : 0.04932 : 0.08636 : 0.04932 : 5 : 7 : 5",
               1);
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
