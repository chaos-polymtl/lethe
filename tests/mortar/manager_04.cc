// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: check generated points using implemented function to
 * read mesh and manifolds, and computing number of faces at the rotor-stator
 * interface and the rotor radius.
 */

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/data_out_faces.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/boundary_conditions.h>
#include <core/grids.h>
#include <core/mortar_coupling_manager.h>
#include <core/parameters.h>

// Tests (with common definitions)
#include <deal.II/base/conditional_ostream.h>

#include <../tests/tests.h>

#include <fstream>

void
test()
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  const unsigned int dim                 = 2;
  const unsigned int mapping_degree      = 3;
  const unsigned int n_quadrature_points = 3;

  Parameters::Mesh                       mesh_parameters;
  Parameters::Mortar<dim>                mortar_parameters;
  Parameters::Manifolds                  manifolds_parameters;
  BoundaryConditions::BoundaryConditions boundary_conditions;
  boundary_conditions.type[0] = BoundaryConditions::BoundaryType::none;

  // Stator mesh parameters
  mesh_parameters.type                     = Parameters::Mesh::Type::dealii;
  mesh_parameters.grid_type                = "hyper_shell";
  mesh_parameters.grid_arguments           = "0, 0 : 0.5 : 1.0 : 6 : true";
  mesh_parameters.scale                    = 1;
  mesh_parameters.simplex                  = false;
  mesh_parameters.initial_refinement       = 2;
  mesh_parameters.refine_until_target_size = false;
  mesh_parameters.boundaries_to_refine     = std::vector<int>();
  mesh_parameters.initial_refinement_at_boundaries = 0;

  // Rotor mesh parameters
  mortar_parameters.enable           = "true";
  mortar_parameters.rotor_mesh       = std::make_shared<Parameters::Mesh>();
  mortar_parameters.rotor_mesh->type = Parameters::Mesh::Type::dealii;
  mortar_parameters.rotor_mesh->grid_type      = "hyper_shell";
  mortar_parameters.rotor_mesh->grid_arguments = "0, 0 : 0.25 : 0.5 : 6 : true";
  mortar_parameters.rotor_mesh->scale          = 1;
  mortar_parameters.rotor_mesh->simplex        = false;
  mortar_parameters.stator_boundary_id         = 0;
  mortar_parameters.rotor_boundary_id          = 3; // after shifting

  // Initialized merged triangulation
  parallel::distributed::Triangulation<dim> triangulation(comm);

  // Merge stator and rotor triangulations
  read_mesh_and_manifolds_for_stator_and_rotor(triangulation,
                                               mesh_parameters,
                                               manifolds_parameters,
                                               false,
                                               boundary_conditions,
                                               mortar_parameters);

  // Number of subdivisions per process
  unsigned int n_subdivisions_local = 0;
  // Number of vertices at the boundary per process
  unsigned int n_vertices_local = 0;
  // Tolerance for rotor radius computation
  const double tolerance = 1e-8;
  // Min and max values for rotor radius computation
  double radius_min = 1e12;
  double radius_max = 1e-12;

  // Check number of faces and vertices at the rotor-stator interface
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (const auto &face : cell->face_iterators())
            {
              if (face->at_boundary())
                {
                  if (face->boundary_id() ==
                      mortar_parameters.rotor_boundary_id)
                    {
                      n_subdivisions_local++;
                      for (unsigned int vertex_index = 0;
                           vertex_index < face->n_vertices();
                           vertex_index++)
                        {
                          n_vertices_local++;
                          auto   v = face->vertex(vertex_index);
                          double radius_current =
                            v.distance(mortar_parameters.center_of_rotation);
                          radius_min = std::min(radius_min, radius_current);
                          radius_max = std::max(radius_max, radius_current);
                        }
                    }
                }
            }
        }
    }

  // Total number of faces
  const unsigned int n_subdivisions =
    Utilities::MPI::sum(n_subdivisions_local, comm);

  // Min and max values over all processes
  radius_min = Utilities::MPI::min(radius_min, comm);
  radius_max = Utilities::MPI::max(radius_max, comm);

  AssertThrow(
    std::abs(radius_max - radius_min) < tolerance,
    ExcMessage(
      "The computed radius of the rotor mesh has a variation greater than the tolerance across the rotor domain, meaning that the prescribed center of rotation and the rotor geometry are not in accordance."));

  // Final radius value
  const double radius = radius_min;

  deallog << "Faces at rotor-stator interface : " << n_subdivisions
          << std::endl;
  deallog << "Rotor radius : " << radius << std::endl;

  // cell angle variation
  const double delta = 2 * numbers::PI / n_subdivisions;

  // rotate inner mesh using random scaling factors
  for (const double scale :
       {0.0, 0.1, 0.5, 0.9, 1.1, 2.1, n_subdivisions - 0.9})
    {
      const double rotate = delta * scale;
      deallog << "Rotation angle: " << rotate << std::endl;

      const MortarManager<dim> manager(n_subdivisions,
                                       n_quadrature_points,
                                       radius,
                                       rotate);

      const auto print = [&](const double shift) {
        for (unsigned int i = 0; i < n_subdivisions; ++i)
          {
            // center point of each cell (in radians)
            const double rad = delta * (i + shift + 0.5);

            deallog << "Shift: " << shift << std::endl;
            deallog << "Cell center at " << rad << ": " << std::endl;

            const auto indices = manager.get_indices(rad);
            const auto weights = manager.get_weights(rad);
            const auto points  = manager.get_points(rad);
            const auto normals = manager.get_normals(rad);

            for (unsigned int i = 0; i < indices.size(); ++i)
              {
                deallog << indices[i] << ": " << weights[i] << "   "
                        << points[i][0] << "   " << points[i][1] << "   "
                        << normals[i][0] << "   " << normals[i][1] << std::endl;
              }
            deallog << std::endl;
          }
        deallog << std::endl;
      };

      // print generated points for aligned mesh (scale = 0.0)
      print(0.0);

      // print generated points for non-aligned mesh
      if (scale != 0.0)
        print(scale - std::floor(scale)); // outer (fixed) mesh

      if (scale != 0.0)
        print(scale); // inner (rotated) mesh

      deallog << std::endl << std::endl << std::endl;
    }

  // Generate vtu file
  DataOut<dim>       data_out;
  MappingQ<dim, dim> mapping(mapping_degree);

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);
  data_out.attach_triangulation(triangulation);

  Vector<double> ranks(triangulation.n_active_cells());
  ranks = Utilities::MPI::this_mpi_process(comm);
  data_out.add_data_vector(ranks, "ranks");
  data_out.build_patches(mapping,
                         mapping_degree + 1,
                         DataOut<dim>::CurvedCellRegion::curved_inner_cells);
  data_out.write_vtu_in_parallel("out.vtu", comm);
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
      test();
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
