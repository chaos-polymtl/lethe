// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <deal.II/base/config.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

using namespace dealii;


/**
 * @brief Create grid with hyper_shell and hyper_shell geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
hyper_shell_with_hyper_shell(const double        radius,
                             Triangulation<dim> &tria,
                             const double        tolerance = 0.0)
{
  const double r1_i = radius * 0.25;
  const double r1_o = radius * 0.5;
  const double r2_i = radius * 0.5;
  const double r2_o = radius * 1.0;

  // inner domain triangulation
  Triangulation<2> circle_one;
  GridGenerator::hyper_shell(circle_one, Point<2>(), r1_i, r1_o, 6, true);
  // outer domain triangulation
  Triangulation<2> circle_two;
  GridGenerator::hyper_shell(circle_two, Point<2>(), r2_i, r2_o, 6, true);

  // shift boundary id of circle two
  for (const auto &face : circle_two.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 2);

  // create unique triangulation
  Triangulation<2> temp;
  GridGenerator::merge_triangulations(
    circle_one, circle_two, temp, tolerance, true, true);
  temp.set_manifold(0, SphericalManifold<2>(Point<2>()));

  if constexpr (dim == 3)
    GridGenerator::extrude_triangulation(temp, 3, radius, tria, true);
  if constexpr (dim == 2)
    tria.copy_triangulation(temp);

  // store manifolds in merged triangulation
  if (dim == 2)
    tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(0, CylindricalManifold<dim>(2));
}

/**
 * @brief Create grid with hyper_ball_balanced and hyper_cube_with_cylindrical_hole geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
hyper_cube_with_cylindrical_hole(const double        radius,
                                 const double        outer_radius,
                                 const double        rotate,
                                 Triangulation<dim> &tria)
{
  Triangulation<2> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

  // outer domain triangulation
  GridGenerator::hyper_cube_with_cylindrical_hole(
    tria_1, radius, outer_radius, 5.0, 1.0, true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  // create unique triangulation
  Triangulation<2> temp;
  GridGenerator::merge_triangulations(tria_0, tria_1, temp, 0, true, true);
  temp.set_manifold(0, SphericalManifold<2>(Point<2>()));
  temp.set_manifold(1, FlatManifold<2>());
  temp.set_manifold(2, SphericalManifold<2>(Point<2>()));

  if constexpr (dim == 3)
    GridGenerator::extrude_triangulation(temp, 2, radius, tria, true);
  if constexpr (dim == 2)
    tria.copy_triangulation(temp);

  // store manifolds in merged triangulation
  if (dim == 2)
    tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(0, CylindricalManifold<dim>(2));
  tria.set_manifold(1, FlatManifold<dim>());
  if (dim == 2)
    tria.set_manifold(2, SphericalManifold<dim>(Point<dim>()));
  else
    tria.set_manifold(2, CylindricalManifold<dim>(2));
}

/**
 * @brief Create grid with hyper_ball_balanced and hyper_cube_with_cylindrical_hole geometries
 * Only the merged triangulation (and resulting grid) are stored
 * Consider a tolerance of 1e-9 when merging boundary points
 */
template <int dim>
void
hyper_cube_with_cylindrical_hole_with_tolerance(const double radius,
                                                const double outer_radius,
                                                const double rotate,
                                                Triangulation<dim> &tria)
{
  Triangulation<dim> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::hyper_ball_balanced(tria_0, {}, radius);
  GridTools::rotate(rotate, tria_0);

  // outer domain triangulation
  GridGenerator::hyper_cube_with_cylindrical_hole(
    tria_1, radius, outer_radius, 5.0, 1.0, true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      {
        face->set_boundary_id(face->boundary_id() + 1);
        face->set_manifold_id(face->manifold_id() + 2);
      }

  // crete unique triangulation
  GridGenerator::merge_triangulations(tria_0, tria_1, tria, 1e-9, true, false);
  // store manifolds in merged triangulation
  tria.set_manifold(0, SphericalManifold<dim>(Point<dim>()));
  tria.set_manifold(1, FlatManifold<dim>());
  tria.set_manifold(2, SphericalManifold<dim>(Point<dim>()));
}

/**
 * @brief Create grid with two subdivided_hyper_rectangle geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
void
split_hyper_cube(Triangulation<2> &tria,
                 const double      left,
                 const double      right,
                 const double      mid,
                 const double      tolerance = 0.0)
{
  Triangulation<2> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::subdivided_hyper_rectangle(tria_0,
                                            std::vector<unsigned int>{1, 1},
                                            Point<2>(left, left),
                                            Point<2>(mid, right),
                                            true);

  // outer domain triangulation
  GridGenerator::subdivided_hyper_rectangle(tria_1,
                                            std::vector<unsigned int>{1, 1},
                                            Point<2>(mid, left),
                                            Point<2>(right, right),
                                            true);

  // shift boundary IDs # in outer grid
  for (const auto &face : tria_1.active_face_iterators())
    if (face->at_boundary())
      face->set_boundary_id(face->boundary_id() + 4);

  // create unique triangulation
  GridGenerator::merge_triangulations(
    tria_0, tria_1, tria, tolerance, true, true);
}

/**
 * @brief Create grid with two subdivided_hyper_rectangle geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
void
split_hyper_cube(Triangulation<1> &tria,
                 const double      left,
                 const double      right,
                 const double      mid,
                 const double      tolerance = 0.0)
{
  Triangulation<1> tria_0, tria_1;

  // inner domain triangulation
  GridGenerator::subdivided_hyper_rectangle(
    tria_0, std::vector<unsigned int>{1}, Point<1>(left), Point<1>(mid), true);

  // outer domain triangulation
  GridGenerator::subdivided_hyper_rectangle(
    tria_1, std::vector<unsigned int>{1}, Point<1>(mid), Point<1>(right), true);

  // shift boundary IDs # in outer grid
  for (const auto &cell : tria_1.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      if (face->at_boundary())
        face->set_boundary_id(face->boundary_id() + 2);

  // create unique triangulation
  GridGenerator::merge_triangulations(
    tria_0, tria_1, tria, tolerance, true, true);
}

/**
 * @brief Create grid with two subdivided_hyper_rectangle geometries
 * Only the merged triangulation (and resulting grid) are stored
 */
template <int dim>
void
split_hyper_cube(Triangulation<dim> &tria,
                 const double        left      = 0.0,
                 const double        right     = 1.0,
                 const double        tolerance = 0.0)
{
  split_hyper_cube(tria, left, right, (left + right) / 2.0, tolerance);
}
