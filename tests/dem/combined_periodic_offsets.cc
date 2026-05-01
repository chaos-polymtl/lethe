// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Unit test for PeriodicBoundariesManipulator::compute_combined_periodic_offsets.
 *
 * Builds a unit hyper_cube triangulation made fully periodic in every
 * direction, drives the public set_periodic_boundaries_information +
 * map_periodic_cells path of PeriodicBoundariesManipulator, then prints
 * the resulting combined_periodic_offsets sorted lexicographically.
 *
 * For a fully periodic unit hyper_cube the offset magnitudes are all 1, so
 * the expected combined_periodic_offsets is exactly the set of non-zero
 * vectors with each component in {-1, 0, +1}: 8 entries in 2D, 26 in 3D.
 *
 * Sorting before printing makes the output independent of insertion order,
 * so the test asserts SET equality (not just size) against the canonical
 * expected set. A wrong-but-same-size list would diff differently.
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

#include <algorithm>

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
  for (unsigned int d = 0; d < dim; ++d)
    GridTools::collect_periodic_faces(triangulation,
                                      /*b_id1=*/2 * d,
                                      /*b_id2=*/2 * d + 1,
                                      /*direction=*/d,
                                      matched_pairs);
  triangulation.add_periodicity(matched_pairs);

  // Two refinements give enough boundary cells with valid periodic
  // neighbours to trigger the offset computation in map_periodic_cells.
  triangulation.refine_global(2);

  // Configure the manipulator with one periodic pair per direction. Use
  // the primary face id (2*d) as both the BC index key and the principal
  // boundary id.
  std::unordered_map<unsigned int, types::boundary_id> primary_ids;
  std::unordered_map<unsigned int, unsigned int>       directions;
  std::vector<unsigned int>                            bc_indices;
  for (unsigned int d = 0; d < dim; ++d)
    {
      primary_ids[d] = 2 * d;
      directions[d]  = d;
      bc_indices.push_back(d);
    }

  PeriodicBoundariesManipulator<dim> manipulator;
  manipulator.set_periodic_boundaries_information(primary_ids,
                                                  directions,
                                                  bc_indices);

  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    cells_info;
  manipulator.map_periodic_cells(triangulation, cells_info);

  // Snapshot and sort lexicographically so the test output is
  // independent of insertion order.
  std::vector<Tensor<1, dim>> offsets =
    manipulator.get_combined_periodic_offsets();
  std::sort(offsets.begin(),
            offsets.end(),
            [](const Tensor<1, dim> &a, const Tensor<1, dim> &b) {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  if (a[d] < b[d])
                    return true;
                  if (a[d] > b[d])
                    return false;
                }
              return false;
            });

  deallog << "dim = " << dim << ", size = " << offsets.size() << std::endl;
  for (const auto &o : offsets)
    {
      deallog << " (";
      for (unsigned int d = 0; d < dim; ++d)
        deallog << (d == 0 ? "" : ", ") << o[d];
      deallog << ")" << std::endl;
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
