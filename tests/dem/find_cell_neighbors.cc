// SPDX-FileCopyrightText: Copyright (c) 2020, 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test generates a triangulation that is twiced refined and check
 * if the cells neighbors are the correct ones.
 */

// Deal.II includes
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <dem/find_cell_neighbors.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim, bool reciprocal>
void
test()
{
  // Creating the mesh and refinement
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  int                                       hyper_cube_length = 1;
  GridGenerator::hyper_cube(triangulation,
                            -1 * hyper_cube_length,
                            hyper_cube_length,
                            true);
  int refinement_number = 2;
  triangulation.refine_global(refinement_number);

  // Finding the cell neighbors
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_local_neighbor_list;
  std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
    cells_ghost_neighbor_list;

  find_cell_neighbors<dim, reciprocal>(triangulation,
                                       cells_local_neighbor_list,
                                       cells_ghost_neighbor_list);

  // Output
  deallog << "reciprocal = "<< reciprocal << std::endl;
  int i = 0;
  for (auto cell = triangulation.begin_active(); cell != triangulation.end();
       ++cell)
    {
      deallog << "neighbors of cell " << cell << " are: ";
      for (auto iterator = cells_local_neighbor_list[i].begin();
           iterator != cells_local_neighbor_list[i].end();
           ++iterator)
        {
          deallog << " " << *iterator;
        }
      deallog << std::endl;

      ++i;
    }
  deallog << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<3, false>();
      test<3, true>();
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
