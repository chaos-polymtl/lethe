// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: check aligned mesh and output points using implemented function to
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
  const unsigned int fe_degree           = 3;
  const unsigned int mapping_degree      = 3;
  const unsigned int n_quadrature_points = 3;
  const unsigned int N                   = 100;

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

  FE_Q<dim>     fe(fe_degree);
  MappingQ<dim> mapping(mapping_degree);
  QGauss<dim>   quadrature(fe_degree + 1);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  // Compute number of subdivisions at the interface and the rotor radius
  const auto [n_subdivisions, radius] =
    compute_n_subdivisions_and_radius(dof_handler, mortar_parameters);

  deallog << "Faces at rotor-stator interface : " << n_subdivisions
          << std::endl;
  deallog << "Rotor radius : " << radius << std::endl;

  // Print information
  for (unsigned int i = 0; i <= N; ++i)
    {
      const double rotate_increment = 2 * numbers::PI / N * i;

      const MortarManager<dim> manager(n_subdivisions,
                                       n_quadrature_points,
                                       radius,
                                       rotate_increment);

      deallog << "Angle : " << rotate_increment << ", Aligned : "
              << static_cast<unsigned int>(manager.is_mesh_aligned())
              << ", # points : " << manager.get_n_points() << std::endl;
    }

  // Generate vtu file
  DataOut<dim> data_out;

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
