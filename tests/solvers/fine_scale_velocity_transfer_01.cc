// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the transfer of the history of the dynamic RBVMS
 * fine-scale velocity across mesh adaptation (FineScaleVelocityTransfer).
 *
 * The history is stored at the Gauss quadrature points of the cells. A
 * tensor-product polynomial of degree k per direction is stored at the
 * (k+1)^dim quadrature points of every cell, for two history entries and dim
 * components. A part of the mesh is then refined and another part is
 * coarsened. Since the history of a refined cell is interpolated to its
 * children with a polynomial of degree k, and the history of coarsened cells
 * is projected to their parent on polynomials of degree k, the transferred
 * values must be equal to the polynomial at the quadrature points of the
 * adapted mesh, up to round-off.
 */

// Deal.II
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>

// Lethe
#include <solvers/fine_scale_velocity_transfer.h>

// Tests
#include <../tests/tests.h>

// C++
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{
  constexpr unsigned int n_history_entries = 2;

  /**
   * @brief Tensor-product polynomial of degree @p degree per direction, which
   * depends on the history entry and on the component.
   */
  template <int dim>
  double
  reference_value(const Point<dim>  &point,
                  const unsigned int degree,
                  const unsigned int history,
                  const unsigned int component)
  {
    double value = 1. + history;
    for (unsigned int d = 0; d < dim; ++d)
      value *= 0.5 * (component + 1) + std::pow(point[d], degree) -
               0.3 * (d + 1) * point[d];
    return value;
  }

  template <int dim>
  void
  check(const unsigned int degree, const unsigned int n_refinements)
  {
    const MPI_Comm     mpi_communicator = MPI_COMM_WORLD;
    const unsigned int n_q_points_1d    = degree + 1;

    parallel::distributed::Triangulation<dim> triangulation(mpi_communicator);
    GridGenerator::hyper_cube(triangulation);
    triangulation.refine_global(n_refinements);

    const QGauss<dim>    quadrature(n_q_points_1d);
    const MappingQ1<dim> mapping;
    const FE_Q<dim>      fe(1);
    FEValues<dim> fe_values(mapping, fe, quadrature, update_quadrature_points);
    const unsigned int n_q_points = quadrature.size();

    // Fill the cell-wise data with the polynomial
    typename FineScaleVelocityTransfer<dim>::CellData cell_data(
      triangulation.n_active_cells());
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          auto &values = cell_data[cell->active_cell_index()];
          for (unsigned int history = 0; history < n_history_entries; ++history)
            for (unsigned int q = 0; q < n_q_points; ++q)
              for (unsigned int d = 0; d < dim; ++d)
                values.emplace_back(reference_value(
                  fe_values.quadrature_point(q), degree, history, d));
        }

    // Refine the cells close to the origin and coarsen those far from it
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          const Point<dim> center = cell->center();
          const double     maximum_coordinate =
            *std::max_element(center.begin_raw(), center.end_raw());
          const double minimum_coordinate =
            *std::min_element(center.begin_raw(), center.end_raw());
          if (maximum_coordinate < 0.5)
            cell->set_refine_flag();
          else if (minimum_coordinate > 0.5)
            cell->set_coarsen_flag();
        }
    triangulation.prepare_coarsening_and_refinement();

    const unsigned int n_cells_before = triangulation.n_global_active_cells();

    FineScaleVelocityTransfer<dim> transfer(triangulation,
                                            mpi_communicator,
                                            n_q_points_1d);
    transfer.prepare_for_coarsening_and_refinement(cell_data);
    triangulation.execute_coarsening_and_refinement();
    const auto transferred_cell_data =
      transfer.unpack_after_coarsening_and_refinement();

    // Compare the transferred values with the polynomial on the adapted mesh
    double       max_error            = 0.;
    unsigned int n_wrong_size_vectors = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          const auto &values = transferred_cell_data[cell->active_cell_index()];
          if (values.size() != n_history_entries * n_q_points * dim)
            {
              ++n_wrong_size_vectors;
              continue;
            }
          unsigned int value = 0;
          for (unsigned int history = 0; history < n_history_entries; ++history)
            for (unsigned int q = 0; q < n_q_points; ++q)
              for (unsigned int d = 0; d < dim; ++d)
                max_error =
                  std::max(max_error,
                           std::abs(
                             values[value++] -
                             reference_value(fe_values.quadrature_point(q),
                                             degree,
                                             history,
                                             d)));
        }
    max_error = Utilities::MPI::max(max_error, mpi_communicator);
    n_wrong_size_vectors =
      Utilities::MPI::sum(n_wrong_size_vectors, mpi_communicator);

    const double tolerance = 1e-12;
    deallog << dim << "D, degree " << degree << std::endl;
    deallog << "  Number of cells before and after adaptation: "
            << n_cells_before << " " << triangulation.n_global_active_cells()
            << std::endl;
    deallog << "  Number of cells with a wrong number of values: "
            << n_wrong_size_vectors << std::endl;
    deallog << "  Transfer is exact: "
            << (max_error < tolerance ? "true" : "false") << std::endl;
  }
} // namespace

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      mpi_initlog();

      check<2>(1, 3);
      check<2>(2, 3);
      check<2>(3, 3);
      check<3>(1, 2);
      check<3>(2, 2);
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
