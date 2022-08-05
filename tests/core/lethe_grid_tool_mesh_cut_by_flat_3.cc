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
* Author: Émile Bergeron , Lucka Barbeau , Polytechnique Montreal, 2021-
*/

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/particles/data_out.h>

// Lethe
#include <core/lethe_grid_tools.h>

#include <deal.II/numerics/data_out.h>

// Tests (with common definitions)
#include <../tests/tests.h>



void
test()
{
  // MPI_Comm mpi_communicator(MPI_COMM_WORLD);


  Triangulation<3> triangulation;

  Triangulation<2, 3>      flat_triangulation;
  std::vector<Point<3>>    vertices_of_flat(GeometryInfo<3>::vertices_per_face);
  std::vector<CellData<2>> flat_cell_data(1);

  DoFHandler<3>    dof_handler;
  DoFHandler<2, 3> flat_dof_handler;

  // Mesh
  GridGenerator::hyper_cube(triangulation, -1, 1);

  // Mesh flat
  Point<3> v0(0, 0, 0);
  Point<3> v1(0.5, 0.5, 0.1);
  Point<3> v2(-0.1, -0.2, 0.7);
  Point<3> v3(0.6, 0.7, 0.8);

  vertices_of_flat[0] = v0;
  vertices_of_flat[1] = v1;
  vertices_of_flat[2] = v2;
  vertices_of_flat[3] = v3;


  flat_cell_data[0].vertices[0] = 0;
  flat_cell_data[0].vertices[1] = 1;
  flat_cell_data[0].vertices[2] = 2;
  flat_cell_data[0].vertices[3] = 3;

  flat_cell_data[0].material_id = 0;
  flat_triangulation.create_triangulation(vertices_of_flat,
                                          flat_cell_data,
                                          SubCellData());
  triangulation.refine_global(3);

  // Attach triangulation to dof_handler

  dof_handler.reinit(triangulation);
  flat_dof_handler.reinit(flat_triangulation);

  Vector<double> subdomain(triangulation.n_active_cells());
  std::map<unsigned int, std::set<typename DoFHandler<3>::active_cell_iterator>>
    vertice_to_cell;

  LetheGridTools::vertices_cell_mapping(dof_handler, vertice_to_cell);
  std::vector<typename DoFHandler<3>::active_cell_iterator> cells_cut;

  const auto &flat_cell = flat_dof_handler.active_cell_iterators().begin();
  cells_cut = LetheGridTools::find_cells_around_flat_cell(dof_handler,
                                                          flat_cell,
                                                          vertice_to_cell);

  std::sort(cells_cut.begin(),
            cells_cut.end(),
            [](typename DoFHandler<3>::active_cell_iterator cell_1,
               typename DoFHandler<3>::active_cell_iterator cell_2) {
              return cell_1->global_active_cell_index() <
                     cell_2->global_active_cell_index();
            });


  for (unsigned int i = 0; i < cells_cut.size(); ++i)
    {
      cells_cut[i]->set_subdomain_id(1);
      subdomain(cells_cut[i]->global_active_cell_index()) =
        cells_cut[i]->global_active_cell_index();
      deallog << "The cell with ID : "
              << cells_cut[i]->global_active_cell_index() << " is cut "
              << std::endl;
    }

    // Printing the final position for all the vertices

#if (DEAL_II_VERSION_MAJOR < 10 && DEAL_II_VERSION_MINOR < 4)
  Legacy::DataOut<3>                   data_out;
  Legacy::DataOut<2, DoFHandler<2, 3>> flat_data_out;
#else
  DataOut<3>    data_out;
  DataOut<2, 3> flat_data_out;
#endif


  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(subdomain, "subdomain");
  data_out.build_patches();
  std::ofstream output("solution.vtu");
  data_out.write_vtu(output);

  flat_data_out.attach_dof_handler(flat_dof_handler);
  flat_data_out.build_patches();
  std::ofstream flat_output("flat_trig.vtu");
  flat_data_out.write_vtu(flat_output);
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
