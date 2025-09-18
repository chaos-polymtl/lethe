// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test generates a triangulation that is refined twice and check
 * if the cells neighbors with repetition are the correct ones.
 */

// Deal.II includes
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

// Lethe
#include <dem/data_containers.h>
#include <dem/find_cell_neighbors.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
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

  // Neighbor cells with repetition
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
    cells_total_neighbor_list;

  find_full_cell_neighbors<dim>(triangulation, cells_total_neighbor_list);

  // Output
  int i = 0;
  for (auto const &[cell_id, cell_neighbors] : cells_total_neighbor_list)
    {
      deallog << "neighbors of cell " << cell_id << " are: ";
      for (auto &iterator : cell_neighbors)
        {
          deallog << " " << iterator->global_active_cell_index();
        }
      deallog << std::endl;

      ++i;
    }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
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
