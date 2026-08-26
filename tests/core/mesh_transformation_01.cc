// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the mesh transformation (scaling, rotation and translation)
 * applied by attach_grid_to_triangulation.
 *
 * The same transformation is applied to the same deal.II grid read as a
 * quad/hex mesh and as a simplex mesh, and to a mesh refinement box built by
 * build_refinement_box_triangulation. All three must place the mesh at the
 * same location, since the transformation is a property of the mesh
 * parameters, not of the type of element or of the use made of the mesh.
 */

// Deal.II
#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <core/grids.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <numbers>

/**
 * @brief Build the parameters of a unit square deal.II mesh which is scaled,
 * rotated by a quarter turn and translated.
 *
 * The unit square is first scaled to [0,2]x[0,2], then rotated
 * counter-clockwise around the origin to [-2,0]x[0,2] and finally translated by
 * one unit along x to [-1,1]x[0,2].
 *
 * @param[in] simplex Whether the mesh is converted to a simplex mesh
 */
Parameters::Mesh
transformed_unit_square_mesh_parameters(const bool simplex)
{
  Parameters::Mesh mesh_parameters;
  mesh_parameters.type           = Parameters::Mesh::Type::dealii;
  mesh_parameters.grid_type      = "hyper_cube";
  mesh_parameters.grid_arguments = "0 : 1 : true";

  mesh_parameters.initial_refinement                  = 1;
  mesh_parameters.initial_refinement_at_boundaries    = 0;
  mesh_parameters.boundaries_to_refine                = {};
  mesh_parameters.refine_until_target_size            = false;
  mesh_parameters.target_size                         = 1.;
  mesh_parameters.simplex                             = simplex;
  mesh_parameters.check_for_diamond_cells             = false;
  mesh_parameters.expand_particle_wall_contact_search = false;

  mesh_parameters.translation    = Tensor<1, 3>({1., 0., 0.});
  mesh_parameters.rotation_axis  = Tensor<1, 3>({0., 0., 1.});
  mesh_parameters.rotation_angle = std::numbers::pi / 2.;
  mesh_parameters.scale          = 2.;

  return mesh_parameters;
}

/**
 * @brief Snap to zero the coordinates which are zero up to round-off, so that
 * the output does not depend on the floating point representation of the sine
 * and the cosine of the rotation angle.
 */
double
clean(const double coordinate)
{
  return std::abs(coordinate) < 1e-12 ? 0. : coordinate;
}

/**
 * @brief Report the number of cells and the bounding box of a triangulation.
 */
void
report_bounding_box(const std::string      &case_name,
                    const Triangulation<2> &triangulation)
{
  const BoundingBox<2> bounding_box =
    GridTools::compute_bounding_box(triangulation);

  const Point<2> bottom_left = bounding_box.get_boundary_points().first;
  const Point<2> top_right   = bounding_box.get_boundary_points().second;

  deallog << case_name << " : " << triangulation.n_global_active_cells()
          << " cells, bounding box [" << clean(bottom_left[0]) << ","
          << clean(top_right[0]) << "]x[" << clean(bottom_left[1]) << ","
          << clean(top_right[1]) << "]" << std::endl;
}

void
test_quad_mesh_transformation()
{
  const Parameters::Mesh mesh_parameters =
    transformed_unit_square_mesh_parameters(false);

  Triangulation<2> triangulation;
  attach_grid_to_triangulation(triangulation, mesh_parameters);
  triangulation.refine_global(mesh_parameters.initial_refinement);

  report_bounding_box("quad mesh", triangulation);
}

void
test_simplex_mesh_transformation()
{
  const Parameters::Mesh mesh_parameters =
    transformed_unit_square_mesh_parameters(true);

  // A simplex mesh is stored in a fully distributed triangulation, since it is
  // built from a description of an already partitioned triangulation. The
  // initial refinement is applied to the quad mesh the simplex mesh is
  // generated from, within attach_grid_to_triangulation.
  parallel::fullydistributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  attach_grid_to_triangulation(triangulation, mesh_parameters);

  report_bounding_box("simplex mesh", triangulation);
}

void
test_refinement_box_transformation()
{
  const Parameters::Mesh box_mesh_parameters =
    transformed_unit_square_mesh_parameters(false);

  Triangulation<2> box_triangulation;
  build_refinement_box_triangulation(box_mesh_parameters, box_triangulation);

  report_bounding_box("refinement box", box_triangulation);
}

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      initlog();

      test_quad_mesh_transformation();
      test_simplex_mesh_transformation();
      test_refinement_box_transformation();
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
