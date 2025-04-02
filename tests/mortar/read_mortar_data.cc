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
 #include <deal.II/grid/grid_generator.h>
 #include <deal.II/grid/grid_tools.h>
 #include <deal.II/grid/manifold_lib.h>
 #include <deal.II/grid/tria.h>
 #include <deal.II/numerics/data_out.h>
 #include <deal.II/fe/mapping_q.h>
 #include <deal.II/numerics/vector_tools.h>
 #include <deal.II/grid/grid_out.h>

// Lethe
#include <core/grids.h>
#include <core/parameters.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <deal.II/base/conditional_ostream.h>
#include <fstream>

void
test()
{
  const MPI_Comm comm = MPI_COMM_WORLD;

  const unsigned int dim                  = 2;
  const unsigned int n_global_refinements = 2;
  const double       radius               = 1.0;
  const double r1_i                       = radius * 0.25;
  const double r1_o                       = radius * 0.5;
  const double r2_i                       = radius * 0.5;
  const double r2_o                       = radius * 1.0;
  const unsigned int mapping_degree       = 3;

  Parameters::Mesh stator_mesh;
  Parameters::Mortar<dim> mortar;
  Parameters::Manifolds manifolds_parameters;
  BoundaryConditions::BoundaryConditions boundary_conditions;
  
  // Stator mesh parameters
  stator_mesh.type = Parameters::Mesh::Type::dealii;
  stator_mesh.grid_type = "hyper_shell";
  stator_mesh.grid_arguments = "0, 0 : 0.5 : 1.0 : 6 : true";
  stator_mesh.scale = 1;
  stator_mesh.simplex = false;
  stator_mesh.initial_refinement = n_global_refinements;
  stator_mesh.refine_until_target_size = false;

  // Rotor mesh parameters
  mortar.enable = "true";
  mortar.rotor_mesh = std::make_shared<Parameters::Mesh>();
  mortar.rotor_mesh->type = Parameters::Mesh::Type::dealii;
  mortar.rotor_mesh->grid_type = "hyper_shell";
  mortar.rotor_mesh->grid_arguments = "0, 0 : 0.25 : 0.5 : 6 : true";
  mortar.rotor_mesh->scale = 1;
  mortar.rotor_mesh->simplex = false;

  parallel::distributed::Triangulation<dim> merged_tria(comm);
  
  // Stator triangulation
  Triangulation<dim> stator_temp_tria;
  attach_grid_to_triangulation(stator_temp_tria, stator_mesh);
  deallog << "Stator boundary IDs :" << stator_temp_tria.get_boundary_ids() << std::endl;
  
  // Rotor triangulation
  Triangulation<dim> rotor_temp_tria;
  attach_grid_to_triangulation(rotor_temp_tria, *mortar.rotor_mesh);
  deallog << "Rotor boundary IDs :" << rotor_temp_tria.get_boundary_ids() << std::endl;

  // Shift rotor boundary IDs #
  for (const auto &face : rotor_temp_tria.active_face_iterators())
    if (face->at_boundary())
    {
      face->set_boundary_id(face->boundary_id() + stator_temp_tria.get_boundary_ids().size());
    }
  deallog << "After Shifting IDs #" << std::endl;
  deallog << "Stator boundary IDs :" << stator_temp_tria.get_boundary_ids() << std::endl;
  deallog << "Rotor boundary IDs :" << rotor_temp_tria.get_boundary_ids() << std::endl;
  
  // Merge triangulations
  GridGenerator::merge_triangulations(stator_temp_tria, 
                                      rotor_temp_tria, 
                                      merged_tria, 
                                      1.0e-12, 
                                      true, 
                                      true);

  // Set manifold                                    
  merged_tria.set_manifold(0,
    SphericalManifold<dim>((dim == 2) ? Point<dim>(0, 0) :
                                        Point<dim>(0, 0, 0)));

  // Initial refinement                                  
  merged_tria.refine_global(n_global_refinements);

  // Print information
  deallog << "Merged triangulation" << std::endl;
  deallog << "Number of active cells : " << merged_tria.n_active_cells() << std::endl;
  deallog << "Number of vertices : " << merged_tria.n_vertices() << std::endl;

  for (const auto &face : merged_tria.active_face_iterators())
    deallog << "Cell center : " << face->center() << std::endl;

  // Generate vtu file
  DataOut<dim> data_out;
  MappingQ<dim, dim> mapping(mapping_degree);

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);
  data_out.attach_triangulation(merged_tria);

  Vector<double> ranks(merged_tria.n_active_cells());
  ranks = Utilities::MPI::this_mpi_process(comm);
  data_out.add_data_vector(ranks, "ranks");
  data_out.build_patches(mapping,
                         mapping_degree + 1,
                         DataOut<dim>::CurvedCellRegion::curved_inner_cells);
  data_out.write_vtu_in_parallel("out.vtu", MPI_COMM_WORLD);
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
