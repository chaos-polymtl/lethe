// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests LetheGridTools::compute_periodic_translations().
 *
 * A distributed rectangle (2D) and box (3D) are made periodic in two
 * directions and refined. The translation of each pair of periodic boundaries
 * is computed from the coarse mesh and printed. The test also checks that
 * every process obtains the same translations, including the processes that
 * do not own any cell at the periodic boundaries.
 */

// Deal.II includes
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>

// Lethe
#include <core/grids.h>
#include <core/lethe_grid_tools.h>
#include <core/periodic_boundary.h>

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

namespace
{
  /**
   * @brief Format a number in fixed notation.
   *
   * @param[in] value Number to format.
   *
   * @return Formatted number.
   */
  std::string
  format(const double value)
  {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << value;
    return stream.str();
  }

  /**
   * @brief Compute and print the periodic translations of a subdivided
   * hyper-rectangle with colorized boundaries.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] repetitions Number of coarse cells in each direction.
   *
   * @param[in] lower_point Lower corner of the hyper-rectangle.
   *
   * @param[in] upper_point Upper corner of the hyper-rectangle.
   *
   * @param[in] periodic_boundaries Pairs of periodic boundaries.
   */
  template <int dim>
  void
  test(const std::vector<unsigned int>      &repetitions,
       const Point<dim>                     &lower_point,
       const Point<dim>                     &upper_point,
       const Parameters::PeriodicBoundaries &periodic_boundaries)
  {
    const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

    parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
    GridGenerator::subdivided_hyper_rectangle(
      triangulation, repetitions, lower_point, upper_point, true);
    setup_periodic_boundary_conditions(triangulation, periodic_boundaries);
    triangulation.refine_global(2);

    const std::vector<LetheGridTools::PeriodicTranslation<dim>> translations =
      LetheGridTools::compute_periodic_translations(triangulation,
                                                    periodic_boundaries);

    deallog << "dim = " << dim
            << ", number of translations: " << translations.size() << std::endl;

    for (const auto &translation : translations)
      {
        deallog << "  boundary id " << translation.boundary_id << ", direction "
                << translation.direction << std::endl;
        deallog << "    offset               : (";
        for (unsigned int d = 0; d < dim; ++d)
          deallog << format(translation.offset[d])
                  << ((d + 1 < dim) ? ", " : ")");
        deallog << std::endl;
        deallog << "    principal coordinate : "
                << format(translation.principal_coordinate) << std::endl;
        deallog << "    neighbor coordinate  : "
                << format(translation.neighbor_coordinate) << std::endl;

        // Every process must obtain exactly the same translation.
        bool identical_on_all_processes = true;
        for (const double value : {translation.offset[translation.direction],
                                   translation.principal_coordinate,
                                   translation.neighbor_coordinate})
          identical_on_all_processes =
            identical_on_all_processes &&
            (Utilities::MPI::min(value, mpi_communicator) ==
             Utilities::MPI::max(value, mpi_communicator));
        deallog << "    identical on all processes : "
                << (identical_on_all_processes ? "true" : "false") << std::endl;
      }
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      // 2D: periodic in x (boundaries 0 and 1) and in y (boundaries 2 and 3)
      {
        Parameters::PeriodicBoundaries periodic_boundaries;
        periodic_boundaries[0] = {1, 0};
        periodic_boundaries[2] = {3, 1};
        test<2>({3, 1},
                Point<2>(-1., 0.5),
                Point<2>(2., 1.5),
                periodic_boundaries);
      }

      // 3D: periodic in x (boundaries 0 and 1) and in z (boundaries 4 and 5)
      {
        Parameters::PeriodicBoundaries periodic_boundaries;
        periodic_boundaries[0] = {1, 0};
        periodic_boundaries[4] = {5, 2};
        test<3>({2, 3, 4},
                Point<3>(0., 0., 0.),
                Point<3>(1., 2., 3.),
                periodic_boundaries);
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
