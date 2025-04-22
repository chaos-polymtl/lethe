// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Mortar: Create points on intersected mesh and determine owners on both sides.
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
 #include <core/mortar_coupling_manager.h>
 
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
 
   const unsigned int dim            = 2;
   const unsigned int mapping_degree = 3;
   const unsigned int n_quadrature_points = 3;
   const double radius = 1.0;
 
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
 
  // Mortar manager
  const MortarManager<dim> mm(4 * Utilities::pow(2, mesh_parameters.initial_refinement + 1),
                              n_quadrature_points,
                              radius,
                              mortar_parameters.rotor_mesh->rotation_angle);

  const unsigned int n_points = mm.get_n_points();

  // convert local/ghost points to indices
  std::vector<double>                  local_values;
  std::vector<types::global_dof_index> is_local;
  std::vector<types::global_dof_index> is_ghost;

  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      for (const auto &face : cell->face_iterators())
        if ((face->boundary_id() == mortar_parameters.rotor_boundary_id) || (face->boundary_id() == mortar_parameters.stator_boundary_id))
          {
            const auto indices = mm.get_indices(point_to_rad(face->center()));
            const auto points  = mm.get_points(point_to_rad(face->center()));
            
            for (unsigned int ii = 0; ii < indices.size(); ++ii)
              {
                unsigned int i = indices[ii];
                unsigned int id_local, id_ghost;

                if (face->boundary_id() == mortar_parameters.rotor_boundary_id)
                  {
                    id_local = i;
                    id_ghost = i + n_points;
                  }
                else if (face->boundary_id() == mortar_parameters.stator_boundary_id)
                  {
                    id_local = i + n_points;
                    id_ghost = i;
                  }

                is_local.emplace_back(id_local);
                is_ghost.emplace_back(id_ghost);

                for (unsigned int d = 0; d < dim; ++d)
                  local_values.emplace_back(points[ii][d]);
              }
          }
  
  // Print information
  for (unsigned int processor_number = 0; processor_number < n_mpi_processes;
    ++processor_number)
 {
   MPI_Barrier(comm);
   if (processor_number == this_mpi_process)
     {
       deallog << "MPI=" << this_mpi_process << std::endl;
       deallog << "# points: " << n_points << std::endl;
       deallog << "Local points: " << std::endl << is_local << std::endl;
       deallog << "Ghost points: " << std::endl << is_ghost << std::endl;
     }
   MPI_Barrier(comm);
 }

  Utilities::MPI::NoncontiguousPartitioner partitioner;
  partitioner.reinit(is_local, is_ghost, comm);

  std::vector<double> ghost_values(local_values.size());
  partitioner.template export_to_ghosted_array<double, dim>(local_values,
                                                            ghost_values);
  
  for (unsigned int i = 0; i < local_values.size(); ++i)
    AssertThrow(std::abs(local_values[i] - ghost_values[i]) < 1e-8,
                ExcInternalError());

  ghost_values.assign(local_values.size(), 0.0);

  partitioner.template export_to_ghosted_array<double, 0>(local_values,
                                                          ghost_values,
                                                          dim);

  for (unsigned int i = 0; i < local_values.size(); ++i)
    AssertThrow(std::abs(local_values[i] - ghost_values[i]) < 1e-8,
                ExcInternalError());
 
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
 