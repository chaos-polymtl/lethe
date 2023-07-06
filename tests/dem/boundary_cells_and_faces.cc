/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Shahab Golshan, Polytechnique Montreal, 2019-
 */

/**
 * @brief This test finds the boundary cells in a system and reports
 * the corresponding information of these cells.
 */

// Deal.II includes
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// Lethe
#include <core/parameters_lagrangian.h>

#include <dem/find_boundary_cells_information.h>

// Tests
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Creating the triangulation and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);
  std::vector<unsigned int> outlet_boundaries;

  // Fining boundary cellds information
  BoundaryCellsInformation<dim> boundary_cells_object;
  boundary_cells_object.build(
    triangulation,
    outlet_boundaries,
    false,
    ConditionalOStream(std::cout,
                       Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0));

  // Reporting the information of boundary cells
  for (auto boundary_cells_information_iterator =
         boundary_cells_object.get_boundary_cells_information().begin();
       boundary_cells_information_iterator !=
       boundary_cells_object.get_boundary_cells_information().end();
       ++boundary_cells_information_iterator)
    {
      auto boundary_information = boundary_cells_information_iterator->second;
      deallog << "Cell " << boundary_information.cell
              << " is on system boundaries (boundary"
              << boundary_information.global_face_id << ")" << std::endl;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test<3>();
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
