// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the Moe scalar limiter on discontinuous Galerkin
 * elements of degree one to three on a mesh made of non-affine cells. It
 * verifies that:
 * - the limiter conserves the integral of an under-resolved solution, modifies
 *   it and does not create a new extremum;
 * - a variation of the solution within a cell that is at the level of round-off
 *   is not amplified when a neighboring cell has a different value.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/numerics/vector_tools.h>

// Lethe
#include <core/vector.h>

#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

/**
 * @brief Under-resolved field made of a sharp front and of oscillations.
 */
template <int dim>
class UnderResolvedField : public Function<dim>
{
public:
  /**
   * @brief Value of the field.
   *
   * @param[in] p Location at which the field is evaluated.
   * @param[in] component Component of the field, unused since it is a scalar.
   *
   * @return Value of the field at the location.
   */
  double
  value(const Point<dim> &p, const unsigned int component = 0) const override
  {
    (void)component;
    return std::tanh((p[0] - 0.5) / 0.02) +
           0.3 * std::sin(40. * p[0]) * std::sin(30. * p[1]);
  }
};

/**
 * @brief Integrate a solution over the domain.
 *
 * @param[in] dof_handler DoFHandler of the solution.
 * @param[in] mapping Mapping of the solution.
 * @param[in] quadrature Quadrature used to integrate.
 * @param[in] solution Solution to integrate.
 *
 * @return Integral of the solution.
 */
template <int dim>
double
integrate_solution(const DoFHandler<dim>  &dof_handler,
                   const Mapping<dim>     &mapping,
                   const Quadrature<dim>  &quadrature,
                   const GlobalVectorType &solution)
{
  FEValues<dim>       fe_values(mapping,
                          dof_handler.get_fe(),
                          quadrature,
                          update_values | update_JxW_values);
  std::vector<double> values(quadrature.size());

  double integral = 0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);
        fe_values.get_function_values(solution, values);
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          integral += values[q] * fe_values.JxW(q);
      }
  return integral;
}

/**
 * @brief Calculate the minimum and the maximum of the locally owned values of
 * a solution.
 *
 * @param[in] solution Solution of which the extrema are calculated.
 * @param[in] locally_owned_dofs Locally owned degrees of freedom.
 *
 * @return Pair made of the minimum and of the maximum of the solution.
 */
std::pair<double, double>
calculate_extrema(const GlobalVectorType &solution,
                  const IndexSet         &locally_owned_dofs)
{
  double min_value = std::numeric_limits<double>::max();
  double max_value = std::numeric_limits<double>::lowest();
  for (const auto i : locally_owned_dofs)
    {
      min_value = std::min(min_value, static_cast<double>(solution(i)));
      max_value = std::max(max_value, static_cast<double>(solution(i)));
    }
  return {min_value, max_value};
}

/**
 * @brief Apply the Moe limiter to two solutions and verify its properties.
 *
 * @param[in] degree Degree of the discontinuous Galerkin element.
 */
template <int dim>
void
test(const unsigned int degree)
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  // Tolerance used for all the comparisons of the test
  const double tolerance = 1e-12;

  // Mesh made of non-affine cells
  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::subdivided_hyper_cube(tria, 8, 0., 1.);
  GridTools::transform(
    [](const Point<dim> &p) {
      Point<dim> q = p;
      for (unsigned int d = 0; d < dim; ++d)
        q[d] += 0.05 * std::sin(2. * numbers::PI * p[(d + 1) % dim]);
      return q;
    },
    tria);

  const FE_DGQ<dim>   fe(degree);
  const MappingQ<dim> mapping(1);
  const QGauss<dim>   quadrature(degree + 1);

  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  const IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs();
  const IndexSet locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(dof_handler);

  GlobalVectorType locally_owned_solution(locally_owned_dofs, mpi_communicator);
  GlobalVectorType locally_relevant_solution(locally_owned_dofs,
                                             locally_relevant_dofs,
                                             mpi_communicator);

  deallog << "Dimension " << dim << " degree " << degree << std::endl;

  // 1. Under-resolved solution
  {
    VectorTools::interpolate(mapping,
                             dof_handler,
                             UnderResolvedField<dim>(),
                             locally_owned_solution);
    locally_relevant_solution = locally_owned_solution;

    const GlobalVectorType initial_solution = locally_owned_solution;
    const auto [initial_min, initial_max] =
      calculate_extrema(initial_solution, locally_owned_dofs);
    const double initial_integral = integrate_solution(
      dof_handler, mapping, quadrature, locally_relevant_solution);

    moe_scalar_limiter<dim>(dof_handler,
                            mapping,
                            quadrature,
                            locally_relevant_solution,
                            locally_owned_solution);

    const double limited_integral = integrate_solution(
      dof_handler, mapping, quadrature, locally_relevant_solution);

    GlobalVectorType difference = locally_owned_solution;
    difference -= initial_solution;
    const auto [limited_min, limited_max] =
      calculate_extrema(locally_owned_solution, locally_owned_dofs);

    deallog << "  Integral is conserved: "
            << (std::abs(limited_integral - initial_integral) < tolerance)
            << std::endl;
    deallog << "  Solution is modified: " << (difference.linfty_norm() > 1e-3)
            << std::endl;
    deallog << "  No new extremum: "
            << (limited_max < initial_max + tolerance &&
                limited_min > initial_min - tolerance)
            << std::endl;
  }

  // 2. Piecewise constant solution with a jump in which the cells of the upper
  // side have a variation of the order of round-off.
  {
    const double lower_value = 0.;
    const double upper_value = 0.5;

    std::vector<types::global_dof_index> local_dof_indices(
      fe.n_dofs_per_cell());
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          const bool is_upper_side = cell->center()[0] > 0.5;
          cell->get_dof_indices(local_dof_indices);
          for (const auto &i : local_dof_indices)
            locally_owned_solution(i) =
              is_upper_side ? upper_value : lower_value;

          // Largest number that is smaller than the upper value
          if (is_upper_side)
            locally_owned_solution(local_dof_indices[0]) =
              std::nextafter(upper_value, lower_value);
        }
    locally_owned_solution.compress(VectorOperation::insert);
    locally_relevant_solution = locally_owned_solution;

    moe_scalar_limiter<dim>(dof_handler,
                            mapping,
                            quadrature,
                            locally_relevant_solution,
                            locally_owned_solution);

    const auto [limited_min, limited_max] =
      calculate_extrema(locally_owned_solution, locally_owned_dofs);

    deallog << "  Round-off variation is not amplified: "
            << (limited_max < upper_value + tolerance &&
                limited_min > lower_value - tolerance)
            << std::endl;
  }
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      for (unsigned int degree = 1; degree <= 3; ++degree)
        test<2>(degree);
      test<3>(2);
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
