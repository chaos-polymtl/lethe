// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test reads two mesh input parameters and merges
 * them in a unique triangulation. The goal is to test the steps on the
 * function read_mesh_and_manifolds only related to reading and merging two
 * distinct triangulations.
 * Based on the test mortar/plot_01.cc
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

  Parameters::Mesh                       mesh_parameters;
  Parameters::Mortar<dim>                mortar_parameters;
  Parameters::Manifolds                  manifolds_parameters;
  BoundaryConditions::BoundaryConditions boundary_conditions;
  boundary_conditions.type[0] = BoundaryConditions::BoundaryType::none;

  // Stator mesh parameters
  mesh_parameters.type                     = Parameters::Mesh::Type::dealii;
  mesh_parameters.grid_type                = "cylinder_shell";
  mesh_parameters.grid_arguments           = "2.0 : 0.5 : 1.0 : 4 : 4 : true";
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
  mortar_parameters.rotor_mesh->grid_type = "cylinder_shell";
  mortar_parameters.rotor_mesh->grid_arguments =
    "2.0 : 0.25 : 0.5 : 4 : 4 : true";
  mortar_parameters.rotor_mesh->scale   = 1;
  mortar_parameters.rotor_mesh->simplex = false;
  mortar_parameters.stator_boundary_id  = 0;
  mortar_parameters.rotor_boundary_id   = 5; // after shifting


  // Initialized merged triangulation
  parallel::distributed::Triangulation<dim> triangulation(comm);

  // Merge stator and rotor triangulations
  read_mesh_and_manifolds_for_stator_and_rotor(triangulation,
                                               mesh_parameters,
                                               manifolds_parameters,
                                               false,
                                               boundary_conditions,
                                               mortar_parameters);

  // Print information
  for (unsigned int processor_number = 0; processor_number < n_mpi_processes;
       ++processor_number)
    {
      MPI_Barrier(comm);
      if (processor_number == this_mpi_process)
        {
          deallog << "MPI=" << this_mpi_process << std::endl;
          deallog << "Number of active cells : "
                  << triangulation.n_active_cells() << std::endl;
          deallog << "Number of vertices : " << triangulation.n_vertices()
                  << std::endl;

          for (const auto &face : triangulation.active_face_iterators())
            deallog << "Cell center : " << face->center() << std::endl;
        }
      MPI_Barrier(comm);
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

  // Plot boundary IDs
  DataPostprocessors::BoundaryIds<dim> boundary_ids;
  DataOutFaces<dim>                    data_out_faces;
  FE_Q<dim>                            dummy_fe(1);

  DoFHandler<dim> dummy_dof_handler(triangulation);
  dummy_dof_handler.distribute_dofs(dummy_fe);

  Vector<double> dummy_solution(dummy_dof_handler.n_dofs());

  data_out_faces.attach_dof_handler(dummy_dof_handler);
  data_out_faces.add_data_vector(dummy_solution, boundary_ids);
  data_out_faces.build_patches();

  std::ofstream out("boundary_ids.vtu");
  data_out_faces.write_vtu(out);
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
