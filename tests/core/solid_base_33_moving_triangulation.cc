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
* Author: Carole-Anne Daunais, Polytechnique Montreal, 2020-
*/

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/data_out.h>

// Lethe
#include <core/parameters.h>
#include <core/solid_base.h>
#include <core/solutions_output.h>
#include <solvers/nitsche.h>

// Tests (with common definitions)
#include <../tests/tests.h>



void
test()
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  SimulationParameters<3> NSparam;
  auto                    param = std::make_shared<Parameters::Nitsche<3>>();
  ParameterHandler        prm;
  std::shared_ptr<parallel::DistributedTriangulationBase<3>> fluid_tria =
    std::make_shared<parallel::distributed::Triangulation<3>>(
      mpi_communicator,
      typename Triangulation<3>::MeshSmoothing(
        Triangulation<3>::smoothing_on_refinement |
        Triangulation<3>::smoothing_on_coarsening));

  std::shared_ptr<parallel::DistributedTriangulationBase<3>> solid_tria =
    std::make_shared<parallel::distributed::Triangulation<3>>(
      mpi_communicator,
      typename Triangulation<3, 3>::MeshSmoothing(
        Triangulation<3, 3>::smoothing_on_refinement |
        Triangulation<3, 3>::smoothing_on_coarsening));

  // Mesh of the solid
  param->solid_mesh.type               = Parameters::Mesh::Type::dealii;
  param->solid_mesh.grid_type          = "hyper_cube";
  param->solid_mesh.grid_arguments     = "-0.5 : 0.5 : false";
  param->solid_mesh.initial_refinement = 3;

  double time_step = 0.01;
  param->solid_velocity.declare_parameters(prm, 3);
  prm.set("Function expression", "-pi*y; pi*x; 0");
  param->solid_velocity.parse_parameters(prm);

  // Mesh of the fluid
  GridGenerator::hyper_cube(*fluid_tria, -1, 1);

  const unsigned int degree_velocity = 1;

  // SolidBase class
  SolidBase<3, 3> solid(param, fluid_tria, degree_velocity);
  solid.initial_setup();
  solid.setup_particles();
  DoFHandler<3, 3> &solid_dh = solid.get_solid_dof_handler();

  DataOut<3> data_out;
  data_out.attach_dof_handler(solid_dh);

  const bool        mapping_all = true;
  const MappingQ<3> mapping(degree_velocity, mapping_all);
  data_out.build_patches(mapping, 1, DataOut<3>::curved_inner_cells);

  PVDHandler pvdhandler;

  data_out.build_patches(mapping, 1, DataOut<3>::curved_inner_cells);
  write_vtu_and_pvd<3>(pvdhandler,
                       data_out,
                       "./",
                       "output_solid_triangulation",
                       0,
                       0,
                       1,
                       mpi_communicator);

  for (unsigned int i = 0; i < 100; ++i)
    {
      solid.move_solid_triangulation(time_step);
      data_out.build_patches(mapping, 1, DataOut<3>::curved_inner_cells);
      double time = (i + 1) * time_step;
      if ((i + 1) % 10 == 0)
        {
          write_vtu_and_pvd<3>(pvdhandler,
                               data_out,
                               "./",
                               "output_solid_triangulation",
                               time,
                               i + 1,
                               1,
                               mpi_communicator);
        }
    }
  // Printing the final position for all the vertices

  const unsigned int n_dofs = solid_dh.n_dofs();
  std::vector<bool>  position_printed(n_dofs, false);

  for (const auto &cell : solid_dh.active_cell_iterators())
    {
      for (unsigned int i = 0; i < GeometryInfo<3>::vertices_per_cell; ++i)
        {
          if (position_printed[cell->vertex_index(i)] == false)
            {
              deallog << "Final position of vertex " << cell->vertex_index(i)
                      << " : " << cell->vertex(i) << std::endl;
              position_printed[cell->vertex_index(i)] = true;
            }
        }
    }
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
