// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: test functions to compute radius, number of
 * subdivisions, and rotate mapping for a 3D case.
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
#include <core/lethe_grid_tools.h>
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
  unsigned int   n_mpi_processes(Utilities::MPI::n_mpi_processes(comm));
  unsigned int   this_mpi_process(Utilities::MPI::this_mpi_process(comm));

  const unsigned int dim            = 3;
  const unsigned int mapping_degree = 3;
  const unsigned int fe_degree      = 3;

  Parameters::Mesh                       mesh_parameters;
  Parameters::Mortar<dim>                mortar_parameters;
  Parameters::Manifolds                  manifolds_parameters;
  BoundaryConditions::BoundaryConditions boundary_conditions;
  boundary_conditions.type[0] = BoundaryConditions::BoundaryType::none;

  // Stator mesh parameters
  mesh_parameters.type                     = Parameters::Mesh::Type::dealii;
  mesh_parameters.grid_type                = "cylinder_shell";
  mesh_parameters.grid_arguments           = "2.0 : 0.5 : 1.0 : 4 : 6 : true";
  mesh_parameters.scale                    = 1;
  mesh_parameters.simplex                  = false;
  mesh_parameters.initial_refinement       = 1;
  mesh_parameters.refine_until_target_size = false;
  mesh_parameters.boundaries_to_refine     = std::vector<int>();
  mesh_parameters.initial_refinement_at_boundaries = 0;
  mesh_parameters.translation                      = Tensor<1, dim>({0, 0, 0});
  mesh_parameters.rotation_axis                    = Tensor<1, dim>({0, 0, 1});
  mesh_parameters.rotation_angle                   = 0.0;

  // Rotor mesh parameters
  mortar_parameters.enable           = "true";
  mortar_parameters.rotor_mesh       = std::make_shared<Parameters::Mesh>();
  mortar_parameters.rotor_mesh->type = Parameters::Mesh::Type::dealii;
  mortar_parameters.rotor_mesh->grid_type = "cylinder_shell";
  mortar_parameters.rotor_mesh->grid_arguments =
    "2.0 : 0.25 : 0.5 : 4 : 6 : true";
  mortar_parameters.rotor_mesh->scale          = 1;
  mortar_parameters.rotor_mesh->simplex        = false;
  mortar_parameters.rotor_mesh->translation    = Tensor<1, dim>({0, 0, 0});
  mortar_parameters.rotor_mesh->rotation_axis  = Tensor<1, dim>({0, 0, 1});
  mortar_parameters.rotor_mesh->rotation_angle = 0.0;
  mortar_parameters.stator_boundary_id         = 0;
  mortar_parameters.rotor_boundary_id          = 5; // after shifting
  mortar_parameters.rotation_axis              = Tensor<1, dim>({0, 0, 1});
  mortar_parameters.center_of_rotation         = Point<dim>();
  const double rotation_angle                  = 0.1;

  // Initialized merged triangulation
  parallel::distributed::Triangulation<dim> triangulation(comm);

  // Merge stator and rotor triangulations
  read_mesh_and_manifolds_for_stator_and_rotor(triangulation,
                                               mesh_parameters,
                                               manifolds_parameters,
                                               false,
                                               boundary_conditions,
                                               mortar_parameters);

  FE_Q<dim>          fe(fe_degree);
  DoFHandler<dim>    dof_handler(triangulation);
  MappingQ<dim, dim> mapping(mapping_degree);
  MappingQCache<dim> mapping_cache(mapping_degree);

  // Distribute dofs
  dof_handler.distribute_dofs(fe);

  // Number of subdivisions and rotor radius
  const auto [n_subdivisions, radius, _] =
    compute_n_subdivisions_and_radius(triangulation,
                                      mapping,
                                      mortar_parameters);

  // Rotate mapping
  LetheGridTools::rotate_mapping(dof_handler,
                                 mapping_cache,
                                 mapping,
                                 radius[0],
                                 rotation_angle,
                                 mortar_parameters.center_of_rotation,
                                 mortar_parameters.rotation_axis);

  // Print information
  if (Utilities::MPI::this_mpi_process(comm) == 0)
    {
      deallog << "Rotation angle (rad) : " << rotation_angle << std::endl;
      deallog << "Radial direction subdivisions: " << n_subdivisions[0]
              << std::endl;
      deallog << "Axial direction subdivisions : " << n_subdivisions[1]
              << std::endl;
      deallog << "Radius : " << radius[0] << std::endl;
      deallog << "Axial direction length : " << radius[1] << std::endl;
    }
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
