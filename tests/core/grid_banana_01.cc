// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// Deal.II
#include <deal.II/base/array_view.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/tensor.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grid_banana.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <array>
#include <cmath>
#include <fstream>
#include <map>

/// Boundary id of the obstacle, as assigned by the colorized grid.
constexpr types::boundary_id banana_boundary_id = 2;


/**
 * @brief Build a banana grid and report the quantities that describe it.
 *
 * @param[in] case_name Name of the case, used for the VTK file it writes.
 *
 * @param[in] grid_arguments Arguments passed to GridBanana.
 *
 * @param[in] global_refinement Number of refinements applied once the grid has
 * been generated, which is what exercises the manifolds it carries.
 */
void
run_test(const std::string &case_name,
         const std::string &grid_arguments,
         const unsigned int global_refinement = 0)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;
  deallog << "Grid arguments: \"" << grid_arguments << "\"" << std::endl;

  Triangulation<3, 3> triangulation;
  GridBanana<3, 3>    grid(grid_arguments);
  grid.make_grid(triangulation);
  triangulation.refine_global(global_refinement);

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;
  deallog << "Number of vertices     : " << triangulation.n_vertices()
          << std::endl;

  deallog << "Mesh volume            : "
          << GridTools::volume(triangulation, MappingQ<3>(2)) << std::endl;

  // A cell folded by the deformation, or by a refinement that the manifolds
  // failed to follow, has a non-positive Jacobian somewhere. Sampling it at the
  // vertices of the cell is not enough, because a trilinear cell can be
  // positive at its eight vertices and folded inside, so the whole reference
  // cell is swept.
  const MappingQ<3>   mapping(1);
  const FE_Nothing<3> finite_element;
  const QIterated<3>  dense_quadrature(QTrapezoid<1>(), 2);
  FEValues<3>         fe_values(mapping,
                        finite_element,
                        dense_quadrature,
                        update_jacobians);

  bool all_cells_are_valid = true;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      fe_values.reinit(cell);
      for (unsigned int q = 0; q < dense_quadrature.size(); ++q)
        all_cells_are_valid =
          all_cells_are_valid &&
          (determinant(Tensor<2, 3>(fe_values.jacobian(q))) > 0.0);
    }

  deallog << "Cells                  : "
          << (all_cells_are_valid ? "all valid" : "SOME ARE FOLDED")
          << std::endl;

  // Count the number of faces per boundary id.
  std::map<types::boundary_id, unsigned int> boundary_face_count;
  for (const auto &cell : triangulation.active_cell_iterators())
    for (const unsigned int f : cell->face_indices())
      if (cell->face(f)->at_boundary())
        boundary_face_count[cell->face(f)->boundary_id()]++;

  for (const auto &[id, count] : boundary_face_count)
    deallog << "Boundary id " << static_cast<int>(id)
            << " face count : " << count << std::endl;

  // Write one VTK file per case for visual checks.
  GridOut           grid_out;
  const std::string vtk_filename = "grid_banana_" + case_name + ".vtk";
  std::ofstream     vtk_out(vtk_filename);
  grid_out.write_vtk(triangulation, vtk_out);
}


/**
 * @brief Build a banana grid, refine it at the obstacle only, and report what
 * that gives.
 *
 * This is the refinement that a simulation asks for through
 * Mesh::initial_refinement_at_boundaries, and the one that leaves hanging nodes
 * behind.
 *
 * @param[in] case_name Name of the case.
 *
 * @param[in] grid_arguments Arguments passed to GridBanana.
 *
 * @param[in] sweeps Number of times the cells touching the obstacle are
 * refined.
 */
void
run_boundary_refinement_test(const std::string &case_name,
                             const std::string &grid_arguments,
                             const unsigned int sweeps)
{
  deallog << "==================================================" << std::endl;
  deallog << "Case: " << case_name << std::endl;

  Triangulation<3, 3> triangulation;
  GridBanana<3, 3>    grid(grid_arguments);
  grid.make_grid(triangulation);

  for (unsigned int sweep = 0; sweep < sweeps; ++sweep)
    {
      for (const auto &cell : triangulation.active_cell_iterators())
        for (const unsigned int f : cell->face_indices())
          if (cell->face(f)->at_boundary() &&
              cell->face(f)->boundary_id() == banana_boundary_id)
            cell->set_refine_flag();
      triangulation.execute_coarsening_and_refinement();
    }

  deallog << "Number of active cells : " << triangulation.n_active_cells()
          << std::endl;

  const MappingQ<3>   mapping(1);
  const FE_Nothing<3> finite_element;
  const QIterated<3>  dense_quadrature(QTrapezoid<1>(), 2);
  FEValues<3>         fe_values(mapping,
                        finite_element,
                        dense_quadrature,
                        update_jacobians);

  bool all_cells_are_valid = true;
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      fe_values.reinit(cell);
      for (unsigned int q = 0; q < dense_quadrature.size(); ++q)
        all_cells_are_valid =
          all_cells_are_valid &&
          (determinant(Tensor<2, 3>(fe_values.jacobian(q))) > 0.0);
    }

  deallog << "Cells                  : "
          << (all_cells_are_valid ? "all valid" : "SOME ARE FOLDED")
          << std::endl;
}


/**
 * @brief Check that the banana manifold maps the ball onto the banana and back.
 *
 * The manifold is only ever asked about points of the obstacle, so the round
 * trip is checked there, over a sphere of the radius of the obstacle.
 *
 * @param[in] case_name Name of the case.
 *
 * @param[in] manifold The manifold to check.
 *
 * @param[in] radius Radius of the sphere the round trip is checked on.
 */
void
check_round_trip(const std::string    &case_name,
                 const BananaManifold &manifold,
                 const double          radius)
{
  double largest_error = 0.0;
  for (unsigned int i = 0; i <= 8; ++i)
    for (unsigned int j = 0; j < 16; ++j)
      {
        const double polar     = numbers::PI * i / 8.0;
        const double azimuthal = 2.0 * numbers::PI * j / 16.0;

        const Point<3> point(radius * std::sin(polar) * std::cos(azimuthal),
                             radius * std::sin(polar) * std::sin(azimuthal),
                             radius * std::cos(polar));

        largest_error =
          std::max(largest_error,
                   (manifold.pull_back(manifold.push_forward(point)) - point)
                     .norm());
      }

  deallog << "Round trip on the obstacle, " << case_name << " : "
          << (largest_error < 1e-12 ? "exact" : "INEXACT") << std::endl;
}


/**
 * @brief Check that the banana manifold keeps new points on the obstacle.
 *
 * This is what the manifold is attached for: a point the manifold builds
 * between two points of the surface of the banana must itself lie on that
 * surface, and not on the straight line between them.
 *
 * @param[in] case_name Name of the case.
 *
 * @param[in] manifold The manifold to check.
 *
 * @param[in] radius Radius of the ball the surface is the image of.
 */
void
check_surface(const std::string    &case_name,
              const BananaManifold &manifold,
              const double          radius)
{
  double largest_error = 0.0;
  for (unsigned int i = 1; i < 8; ++i)
    for (unsigned int j = 0; j < 16; ++j)
      {
        const double polar     = numbers::PI * i / 8.0;
        const double azimuthal = 2.0 * numbers::PI * j / 16.0;

        // Two neighbouring points of the surface, and the point the manifold
        // builds halfway between them.
        const Point<3> first(radius * std::sin(polar) * std::cos(azimuthal),
                             radius * std::sin(polar) * std::sin(azimuthal),
                             radius * std::cos(polar));
        const Point<3> second(
          radius * std::sin(polar + numbers::PI / 8.0) * std::cos(azimuthal),
          radius * std::sin(polar + numbers::PI / 8.0) * std::sin(azimuthal),
          radius * std::cos(polar + numbers::PI / 8.0));

        const std::array<Point<3>, 2> deformed = {
          {manifold.push_forward(first), manifold.push_forward(second)}};
        const std::array<double, 2> weights = {{0.5, 0.5}};

        const Point<3> middle =
          manifold.get_new_point(make_array_view(deformed),
                                 make_array_view(weights));

        largest_error = std::max(largest_error,
                                 std::abs(manifold.pull_back(middle).norm() -
                                          radius));
      }

  deallog << "New points on the obstacle, " << case_name << " : "
          << (largest_error < 1e-10 ? "on the surface" : "OFF the surface")
          << std::endl;
}


int
main()
{
  try
    {
      initlog();

      // A straight banana, which is a prolate spheroid, and the same grid bent
      // into a banana. The two share every other argument, so the difference
      // between them is the bend alone.
      run_test("straight",
               "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0,0 : 0.25 : 1 : 2.6 : true : true");

      run_test("bent",
               "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0,0.4 : 0.25 : 1 : 2.6 : true : true");

      // The same banana, refined once while the obstacle is still a ball,
      // and then the same banana refined once after it has been deformed. The
      // manifolds of the deformed grid describe both its surface and the region
      // around it, so the two must agree: refining the grid generated with a
      // refinement of 1 has to give the grid that a refinement of 2 generates.
      run_test("bent_fine",
               "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0,0.4 : 0.25 : 2 : 2.6 : true : true");

      run_test("bent_refined",
               "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0,0.4 : 0.25 : 1 : 2.6 : true : true",
               1);

      // Refining the cells that touch the obstacle, which leaves hanging
      // nodes, is what a simulation does through initial boundary refinement.
      run_boundary_refinement_test(
        "bent_refined_at_the_obstacle",
        "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0,0.4 : 0.25 : 2 : 2.6 : true : true",
        2);

      // Bending along the streamwise direction instead of across it.
      run_test("bent_streamwise",
               "2,3,2,2,2,2 : 0.5 : 0.6 : 4 : 0.4,0 : 0.25 : 1 : 2.6 : true : true");

      check_round_trip("straight",
                       BananaManifold(4.0, 0.5, Tensor<1, 2>({0.0, 0.0}), 0.25),
                       0.5);
      check_round_trip("bent",
                       BananaManifold(4.0, 0.5, Tensor<1, 2>({0.0, 0.4}), 0.25),
                       0.5);
      check_round_trip("bent, on the shell around the obstacle",
                       BananaManifold(4.0, 0.5, Tensor<1, 2>({0.0, 0.4}), 0.25),
                       0.55);

      check_surface("straight",
                    BananaManifold(4.0, 0.5, Tensor<1, 2>({0.0, 0.0}), 0.25),
                    0.5);
      check_surface("bent",
                    BananaManifold(4.0, 0.5, Tensor<1, 2>({0.0, 0.4}), 0.25),
                    0.5);
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
