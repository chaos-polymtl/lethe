/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 * This test generates a flat plane using a simplex mesh out of the deal.II grid
 * generator and calculates the cells that it cuts on an hyper_ball.
 * The list of cut cells/triangles combo is printed and a vtu file with the
 * intersections is generated.
 */

#include <core/parameters.h>
#include <core/serial_solid.h>
#include <core/solid_objects_parameters.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/data_out.h>

#include <../tests/tests.h>

void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Generate a background fluid triangulation made of a sphere
  std::shared_ptr<parallel::TriangulationBase<3>> background_tria =
    std::make_shared<parallel::distributed::Triangulation<3>>(
      mpi_communicator,
      typename Triangulation<3>::MeshSmoothing(
        Triangulation<3>::smoothing_on_refinement |
        Triangulation<3>::smoothing_on_coarsening));

  // Mesh of the fluid
  GridGenerator::hyper_ball(*background_tria, {0.2, 0, 0}, 2);
  background_tria->refine_global(3);

  // Generate the serial solid

  // Parameters for the Serial solid object
  auto param             = std::make_shared<Parameters::RigidSolidObject<3>>();
  param->solid_mesh.type = Parameters::Mesh::Type::dealii;
  param->solid_mesh.grid_type          = "hyper_rectangle";
  param->solid_mesh.grid_arguments     = "-2, -1 : 2, 1 : false";
  param->solid_mesh.initial_refinement = 3;
  param->solid_mesh.simplex            = true;
  param->solid_mesh.translation        = Tensor<1, 3>({0., 0., 0.});
  param->solid_mesh.rotation_axis      = Tensor<1, 3>({1., 0., 0.});
  param->solid_mesh.rotation_angle     = 0.;

  SerialSolid<2, 3> solid(param, 0);

  // Calculate the intersections between the background triangulation and the
  // floating solid
  auto cell_and_triangle_intersection =
    solid.map_solid_in_background_triangulation(*background_tria);
  deallog << "Cell pairs (background, solid)" << std::endl;
  for (auto &cell_pair : cell_and_triangle_intersection)
    {
      deallog << cell_pair.first->id() << " " << cell_pair.second->id()
              << std::endl;
    }

  // Generate sufficient information to generate a graphical output
  DoFHandler<3> background_dof_handler;
  background_dof_handler.reinit(*background_tria);

  DoFHandler<2, 3> &solid_dof_handler = solid.get_dof_handler();
  std::shared_ptr<Triangulation<2, 3>> solid_triangulation =
    solid.get_triangulation();

  // Loop over background
  Vector<double> subdomain_background(background_tria->n_active_cells());
  Vector<double> subdomain_solid(solid_triangulation->n_active_cells());
  for (auto &cell_pair : cell_and_triangle_intersection)
    {
      cell_pair.first->set_subdomain_id(1);
      cell_pair.second->set_subdomain_id(1);
      subdomain_background(cell_pair.first->global_active_cell_index()) = 1;
      subdomain_solid(cell_pair.second->global_active_cell_index())     = 1;
    }


  DataOut<3>    background_data_out;
  DataOut<2, 3> solid_data_out;

  // Generate a VTU for debugging purposes which shows the intersection
  background_data_out.attach_dof_handler(background_dof_handler);
  background_data_out.add_data_vector(subdomain_background, "subdomain");
  background_data_out.build_patches();
  std::ofstream output("background.vtu");
  background_data_out.write_vtu(output);

  solid_data_out.attach_dof_handler(solid_dof_handler);
  solid_data_out.add_data_vector(subdomain_solid, "subdomain");
  solid_data_out.build_patches();
  std::ofstream solid_output("solid_tria.vtu");
  solid_data_out.write_vtu(solid_output);
}

int
main(int argc, char *argv[])
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
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
