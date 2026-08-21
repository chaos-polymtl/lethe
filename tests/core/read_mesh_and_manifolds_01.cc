// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests read_mesh_and_manifolds, the single mesh reading routine shared
 * by the CFD and the DEM solvers.
 *
 * Four behaviours are exercised:
 *  - the plain initial refinement, and the fact that it is skipped on restart,
 *    since the triangulation of a restarted simulation is subsequently replaced
 *    by the one stored in the checkpoint;
 *  - the refinement until a target cell size, which must be skipped on restart
 *    for the same reason;
 *  - the refinement in the vicinity of a list of boundaries
 *    ("initial boundary refinement" / "boundaries refined");
 *  - the attachment of the manifolds declared in the parameter file.
 */

// Deal.II
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_tools.h>

// Lethe
#include <core/grids.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <algorithm>

/**
 * @brief Build the mesh parameters of a unit square discretized by a single
 * deal.II hyper_cube cell.
 */
Parameters::Mesh
unit_square_mesh_parameters()
{
  Parameters::Mesh mesh_parameters;
  mesh_parameters.type           = Parameters::Mesh::Type::dealii;
  mesh_parameters.grid_type      = "hyper_cube";
  mesh_parameters.grid_arguments = "0 : 1 : true";

  mesh_parameters.initial_refinement                  = 0;
  mesh_parameters.initial_refinement_at_boundaries    = 0;
  mesh_parameters.boundaries_to_refine                = {};
  mesh_parameters.refine_until_target_size            = false;
  mesh_parameters.target_size                         = 1.;
  mesh_parameters.simplex                             = false;
  mesh_parameters.check_for_diamond_cells             = false;
  mesh_parameters.expand_particle_wall_contact_search = false;

  mesh_parameters.translation    = Tensor<1, 3>();
  mesh_parameters.rotation_axis  = Tensor<1, 3>({1., 0., 0.});
  mesh_parameters.rotation_angle = 0.;
  mesh_parameters.scale          = 1.;

  return mesh_parameters;
}

/**
 * @brief Read a mesh with the given parameters and report the resulting number
 * of active cells.
 */
void
report_number_of_cells(
  const std::string           &case_name,
  const Parameters::Mesh      &mesh_parameters,
  const bool                   restart,
  const Parameters::Manifolds &manifolds_parameters = Parameters::Manifolds())
{
  parallel::distributed::Triangulation<2> triangulation(
    MPI_COMM_WORLD, Triangulation<2>::limit_level_difference_at_vertices);

  // No periodic boundary is set up in this test
  const Parameters::PeriodicBoundaryInformation periodic_boundary_information;

  read_mesh_and_manifolds(triangulation,
                          mesh_parameters,
                          manifolds_parameters,
                          restart,
                          periodic_boundary_information);

  deallog << case_name << " : " << triangulation.n_global_active_cells()
          << " cells" << std::endl;
}

void
test_initial_refinement()
{
  deallog << "--- Initial refinement ---" << std::endl;

  Parameters::Mesh mesh_parameters   = unit_square_mesh_parameters();
  mesh_parameters.initial_refinement = 2;

  // Without restart, the mesh is refined twice: 4^2 = 16 cells
  report_number_of_cells("initial refinement, no restart",
                         mesh_parameters,
                         false);

  // On restart, no refinement is applied: the mesh stored in the checkpoint
  // already carries its refinement
  report_number_of_cells("initial refinement, restart", mesh_parameters, true);
}

void
test_refinement_until_target_size()
{
  deallog << "--- Refinement until target size ---" << std::endl;

  Parameters::Mesh mesh_parameters         = unit_square_mesh_parameters();
  mesh_parameters.refine_until_target_size = true;

  // The minimal cell diameter of the single-cell unit square is its diagonal,
  // sqrt(2). The number of refinements is floor(log2(sqrt(2)/target_size)), so
  // a target size of 0.1 gives floor(log2(14.14)) = 3 refinements, i.e.
  // 8^2 = 64 cells.
  mesh_parameters.target_size = 0.1;

  report_number_of_cells("target size, no restart", mesh_parameters, false);

  // On restart the target size must be ignored, otherwise the mesh loaded from
  // the checkpoint would be refined a second time
  report_number_of_cells("target size, restart", mesh_parameters, true);
}

void
test_refinement_at_boundaries()
{
  deallog << "--- Refinement at boundaries ---" << std::endl;

  Parameters::Mesh mesh_parameters   = unit_square_mesh_parameters();
  mesh_parameters.initial_refinement = 2;

  // The hyper_cube is colorized: boundary id 0 is the x = 0 face. Refining it
  // once more adds 3 cells for each of the 4 cells touching that boundary.
  mesh_parameters.boundaries_to_refine             = {0};
  mesh_parameters.initial_refinement_at_boundaries = 1;

  report_number_of_cells("one boundary refined once", mesh_parameters, false);

  // Refining two opposite boundaries twice
  mesh_parameters.boundaries_to_refine             = {0, 1};
  mesh_parameters.initial_refinement_at_boundaries = 2;

  report_number_of_cells("two boundaries refined twice",
                         mesh_parameters,
                         false);

  // The boundary refinement is skipped on restart, like the initial refinement
  report_number_of_cells("two boundaries refined twice, restart",
                         mesh_parameters,
                         true);
}

void
test_manifold_attachment()
{
  deallog << "--- Manifold attachment ---" << std::endl;

  const Parameters::Mesh mesh_parameters = unit_square_mesh_parameters();

  // Declare a spherical manifold centered at the origin and attach it to
  // manifold id 1
  Parameters::Manifolds manifolds_parameters;
  manifolds_parameters.size  = 1;
  manifolds_parameters.id    = {1};
  manifolds_parameters.types = {Parameters::Manifolds::ManifoldType::spherical};
  manifolds_parameters.manifold_point     = {"0,0,0"};
  manifolds_parameters.manifold_direction = {"0,1,0"};
  manifolds_parameters.cad_files          = {""};

  parallel::distributed::Triangulation<2> triangulation(
    MPI_COMM_WORLD, Triangulation<2>::limit_level_difference_at_vertices);

  const Parameters::PeriodicBoundaryInformation periodic_boundary_information;

  read_mesh_and_manifolds(triangulation,
                          mesh_parameters,
                          manifolds_parameters,
                          false,
                          periodic_boundary_information);

  // Probe the manifold attached to id 1 by asking for the point halfway between
  // two points lying on the unit circle. A spherical manifold places it back on
  // the circle (norm 1), whereas the default flat manifold would return their
  // arithmetic mean (norm sqrt(2)/2 = 0.707).
  const Point<2> first_point(1., 0.);
  const Point<2> second_point(0., 1.);

  const Point<2> new_point =
    triangulation.get_manifold(1).get_intermediate_point(first_point,
                                                         second_point,
                                                         0.5);

  deallog << "Manifold id 1 midpoint norm : " << new_point.norm() << std::endl;

  // Reading the same mesh without declaring any manifold must leave the
  // triangulation without a manifold of id 1
  parallel::distributed::Triangulation<2> bare_triangulation(
    MPI_COMM_WORLD, Triangulation<2>::limit_level_difference_at_vertices);

  read_mesh_and_manifolds(bare_triangulation,
                          mesh_parameters,
                          Parameters::Manifolds(),
                          false,
                          periodic_boundary_information);

  const std::vector<types::manifold_id> bare_manifold_ids =
    bare_triangulation.get_manifold_ids();
  const bool manifold_one_is_attached =
    std::ranges::find(bare_manifold_ids, 1) != bare_manifold_ids.end();

  deallog << "Manifold id 1 attached without a manifolds subsection : "
          << (manifold_one_is_attached ? "true" : "false") << std::endl;
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();

      test_initial_refinement();
      test_refinement_until_target_size();
      test_refinement_at_boundaries();
      test_manifold_attachment();
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
