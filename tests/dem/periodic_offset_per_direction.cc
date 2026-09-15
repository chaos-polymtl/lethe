// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Unit test for
 * PeriodicBoundariesManipulator::compute_periodic_offset_per_direction.
 *
 * Builds a unit hyper_cube triangulation made fully periodic in every
 * direction, drives the public set_periodic_boundaries_information +
 * map_periodic_cells path of PeriodicBoundariesManipulator, then prints the
 * resulting periodic_offset_per_direction.
 *
 * For a fully periodic unit hyper_cube the period along every direction is
 * 1, so the expected periodic_offset_per_direction is (1, 1) in 2D and
 * (1, 1, 1) in 3D.
 */

// Deal.II
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// Lethe
#include <dem/data_containers.h>
#include <dem/periodic_boundaries_manipulator.h>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

template <int dim>
void
test()
{
  // Fully periodic unit hyper_cube; colorize so opposite faces have IDs
  // (2*d, 2*d+1) per direction d.
  parallel::distributed::Triangulation<dim> triangulation(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(triangulation, 0., 1., /*colorize=*/true);

  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<dim>::cell_iterator>>
    matched_pairs;
  for (int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(triangulation,
                                      /*b_id1=*/2 * d,
                                      /*b_id2=*/2 * d + 1,
                                      /*direction=*/d,
                                      matched_pairs);
  triangulation.add_periodicity(matched_pairs);

  // Two refinements give enough boundary cells with valid periodic
  // neighbours to trigger the offset computation in map_periodic_cells.
  triangulation.refine_global(2);

  // Configure the manipulator with one periodic pair per direction. The map is
  // keyed by the principal periodic boundary id (2*d) and matches the pairs
  // collected above.
  Parameters::PeriodicBoundaries periodic_boundaries;
  for (int d = 0; d < dim; ++d)
    periodic_boundaries[2 * d] = {.neighbor_id =
                                    static_cast<types::boundary_id>(2 * d + 1),
                                  .direction = static_cast<unsigned int>(d)};

  PeriodicBoundariesManipulator<dim> manipulator;
  manipulator.set_periodic_boundaries_information(periodic_boundaries);

  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    cells_info;
  manipulator.map_periodic_cells(triangulation, cells_info);

  const Tensor<1, dim> &offset_per_direction =
    manipulator.get_periodic_offset_per_direction();

  deallog << "dim = " << dim << ", periodic_offset_per_direction = (";
  for (int d = 0; d < dim; ++d)
    deallog << (d == 0 ? "" : ", ") << offset_per_direction[d];
  deallog << ")" << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      test<2>();
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
