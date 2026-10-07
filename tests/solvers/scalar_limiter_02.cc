// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the Kuzmin scalar limiter on discontinuous Galerkin
 * elements of degree one to three on a mesh made of non-affine cells. It
 * verifies that:
 * - the limiter conserves the integral of an under-resolved solution and
 *   modifies it;
 * - for a degree of one, the limited solution at every vertex that is not on a
 *   boundary is between the minimum and the maximum of the means of the cells
 *   that share this vertex;
 * - for a degree of two and above, a smooth solution with an extremum within
 *   the domain is not modified.
 */

// Deal.II includes
#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
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
 * @brief Polynomial of degree two with a maximum within the domain.
 */
template <int dim>
class SmoothExtremumField : public Function<dim>
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
    double value = 1.;
    for (unsigned int d = 0; d < dim; ++d)
      value -= (p[d] - 0.52) * (p[d] - 0.52);
    return value;
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
 * @brief Verify that the solution at every vertex that is not on a boundary is
 * between the minimum and the maximum of the means of the cells that share
 * this vertex.
 *
 * @param[in] dof_handler DoFHandler of the solution.
 * @param[in] mapping Mapping of the solution.
 * @param[in] quadrature Quadrature used to calculate the means.
 * @param[in] solution Solution to verify.
 * @param[in] tolerance Tolerance on the bounds.
 *
 * @return True if the solution is within the bounds at every vertex.
 */
template <int dim>
bool
vertex_values_are_bounded(const DoFHandler<dim>  &dof_handler,
                          const Mapping<dim>     &mapping,
                          const Quadrature<dim>  &quadrature,
                          const GlobalVectorType &solution,
                          const double            tolerance)
{
  const auto        &triangulation = dof_handler.get_triangulation();
  const unsigned int n_vertices    = triangulation.n_vertices();

  // Mean of every cell and bounds of every vertex
  std::vector<double> mean_per_cell(triangulation.n_active_cells());
  std::vector<double> min_mean_per_vertex(n_vertices,
                                          std::numeric_limits<double>::max());
  std::vector<double> max_mean_per_vertex(
    n_vertices, std::numeric_limits<double>::lowest());
  std::vector<bool> vertex_is_at_boundary(n_vertices, false);

  FEValues<dim>       fe_values(mapping,
                          dof_handler.get_fe(),
                          quadrature,
                          update_values | update_JxW_values);
  std::vector<double> values(quadrature.size());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      fe_values.get_function_values(solution, values);
      double integral = 0;
      double measure  = 0;
      for (unsigned int q = 0; q < quadrature.size(); ++q)
        {
          integral += values[q] * fe_values.JxW(q);
          measure += fe_values.JxW(q);
        }
      const double mean                        = integral / measure;
      mean_per_cell[cell->active_cell_index()] = mean;

      for (unsigned int v = 0; v < cell->n_vertices(); ++v)
        {
          const unsigned int vertex_index = cell->vertex_index(v);
          min_mean_per_vertex[vertex_index] =
            std::min(min_mean_per_vertex[vertex_index], mean);
          max_mean_per_vertex[vertex_index] =
            std::max(max_mean_per_vertex[vertex_index], mean);
        }

      for (const auto f : cell->face_indices())
        if (cell->at_boundary(f))
          for (unsigned int v = 0; v < cell->face(f)->n_vertices(); ++v)
            vertex_is_at_boundary[cell->face(f)->vertex_index(v)] = true;
    }

  // Solution at the vertices of every cell
  std::vector<Point<dim>> unit_vertices;
  for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    unit_vertices.push_back(GeometryInfo<dim>::unit_cell_vertex(v));
  const Quadrature<dim> vertex_quadrature(unit_vertices);
  FEValues<dim>         fe_values_vertices(mapping,
                                   dof_handler.get_fe(),
                                   vertex_quadrature,
                                   update_values);
  std::vector<double>   vertex_values(vertex_quadrature.size());

  bool is_bounded = true;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values_vertices.reinit(cell);
      fe_values_vertices.get_function_values(solution, vertex_values);
      for (unsigned int v = 0; v < cell->n_vertices(); ++v)
        {
          const unsigned int vertex_index = cell->vertex_index(v);
          if (vertex_is_at_boundary[vertex_index])
            continue;

          is_bounded =
            is_bounded &&
            vertex_values[v] < max_mean_per_vertex[vertex_index] + tolerance &&
            vertex_values[v] > min_mean_per_vertex[vertex_index] - tolerance;
        }
    }
  return is_bounded;
}

/**
 * @brief Apply the Kuzmin limiter to two solutions and verify its properties.
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
    const double           initial_integral = integrate_solution(
      dof_handler, mapping, quadrature, locally_relevant_solution);

    kuzmin_scalar_limiter<dim>(dof_handler,
                               mapping,
                               quadrature,
                               locally_relevant_solution,
                               locally_owned_solution);

    const double limited_integral = integrate_solution(
      dof_handler, mapping, quadrature, locally_relevant_solution);

    GlobalVectorType difference = locally_owned_solution;
    difference -= initial_solution;

    deallog << "  Integral is conserved: "
            << (std::abs(limited_integral - initial_integral) < tolerance)
            << std::endl;
    deallog << "  Solution is modified: " << (difference.linfty_norm() > 1e-3)
            << std::endl;

    if (degree == 1)
      deallog << "  Vertex values are within the vertex bounds: "
              << vertex_values_are_bounded(dof_handler,
                                           mapping,
                                           quadrature,
                                           locally_relevant_solution,
                                           tolerance)
              << std::endl;
  }

  // 2. Smooth solution with an extremum. The hierarchical limiting requires a
  // degree of at least two to preserve it.
  if (degree >= 2)
    {
      VectorTools::interpolate(mapping,
                               dof_handler,
                               SmoothExtremumField<dim>(),
                               locally_owned_solution);
      locally_relevant_solution = locally_owned_solution;

      const GlobalVectorType initial_solution = locally_owned_solution;

      kuzmin_scalar_limiter<dim>(dof_handler,
                                 mapping,
                                 quadrature,
                                 locally_relevant_solution,
                                 locally_owned_solution);

      GlobalVectorType difference = locally_owned_solution;
      difference -= initial_solution;

      deallog << "  Smooth extremum is preserved: "
              << (difference.linfty_norm() < tolerance) << std::endl;
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
