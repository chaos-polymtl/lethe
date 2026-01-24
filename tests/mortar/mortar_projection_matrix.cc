// SPDX-FileCopyrightText: Copyright (c) 2025-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Compare interpolation and projection operations onto a subcell
 * of size factor.
 */

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

using namespace dealii;

template <int dim>
Quadrature<dim>
scale(const Quadrature<dim> &quadrature_2, const double factor)
{
  std::vector<Point<dim>> points;
  for (const auto p : quadrature_2.get_points())
    points.emplace_back(p * factor);
  return Quadrature<dim>(points);
}

template <int dim, typename number, int spacedim>
void
get_projection_matrix(const FiniteElement<dim, spacedim> &fe1,
                      const double                        factor,
                      FullMatrix<number>                 &matrix)
{
  matrix = 0;

  const unsigned int n1 = fe1.n_dofs_per_cell();
  const unsigned int n2 = fe1.n_dofs_per_cell();
  const unsigned int nd = fe1.n_components();

  const ReferenceCell reference_cell = fe1.reference_cell();

  // First, create a local mass matrix for the unit cell
  Triangulation<dim, spacedim> tr_1;
  GridGenerator::reference_cell(tr_1, reference_cell);
  Triangulation<dim, spacedim> tr_2;
  GridGenerator::reference_cell(tr_2, reference_cell);
  GridTools::scale(factor, tr_2);

  const auto &mapping =
    reference_cell.template get_default_linear_mapping<dim, spacedim>();

  // Choose a Gauss quadrature rule that is exact up to degree 2n-1
  const unsigned int degree = fe1.tensor_degree();
  Assert(degree != numbers::invalid_unsigned_int, ExcNotImplemented());
  const auto quadrature_2 =
    reference_cell.get_gauss_type_quadrature<dim>(degree + 1);

  std::vector<Point<dim>> points;
  for (const auto p : quadrature_2.get_points())
    points.emplace_back(p * factor);
  Quadrature<dim> quadrature_1(points);

  const unsigned int nq = quadrature_2.size();

  // Set up FEValues.
  const UpdateFlags flags =
    update_values | update_quadrature_points | update_JxW_values;
  FEValues<dim> val1(mapping, fe1, quadrature_1, update_values);
  val1.reinit(tr_1.begin_active());

  FEValues<dim> val2(mapping, fe1, quadrature_2, flags);
  val2.reinit(tr_2.begin_active());

  // Integrate and invert mass matrix. This happens in the target space
  FullMatrix<double> mass(n2, n2);

  for (unsigned int i = 0; i < n2; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      for (unsigned int d = 0; d < nd; ++d)
        for (unsigned int k = 0; k < nq; ++k)
          mass(i, j) += val2.JxW(k) * val2.shape_value_component(i, k, d) *
                        val2.shape_value_component(j, k, d);

  // Invert the matrix. Gauss-Jordan should be sufficient since we expect
  // the mass matrix to be well-conditioned
  mass.gauss_jordan();

  // Now, test every function of fe1 with test functions of fe2 and
  // compute the projection of each unit vector.
  Vector<double> b(n2);
  Vector<double> x(n2);

  for (unsigned int j = 0; j < n1; ++j)
    {
      b = 0.;
      for (unsigned int i = 0; i < n2; ++i)
        for (unsigned int k = 0; k < quadrature_2.size(); ++k)
          for (unsigned int d = 0; d < nd; ++d)
            b(i) += val1.shape_value_component(j, k, d) *
                    val2.shape_value_component(i, k, d) * val2.JxW(k);

      // Multiply by the inverse
      mass.vmult(x, b);
      for (unsigned int i = 0; i < n2; ++i)
        matrix(i, j) = x(i);
    }
}

template <int dim, typename number, int spacedim>
void
get_interpolation_matrix(const FiniteElement<dim, spacedim> &fe1,
                         const Quadrature<dim>              &quadrature,
                         FullMatrix<number>                 &matrix)
{
  matrix = 0;

  const unsigned int n2 = fe1.n_dofs_per_cell();
  const unsigned int nd = fe1.n_components();

  const ReferenceCell reference_cell = fe1.reference_cell();

  // First, create a local mass matrix for the unit cell
  Triangulation<dim, spacedim> tr_2;
  GridGenerator::reference_cell(tr_2, reference_cell);

  const auto &mapping =
    reference_cell.template get_default_linear_mapping<dim, spacedim>();

  const unsigned int nq = quadrature.size();

  // Set up FEValues.
  const UpdateFlags flags =
    update_values | update_quadrature_points | update_JxW_values;

  FEValues<dim> val2(mapping, fe1, quadrature, flags);
  val2.reinit(tr_2.begin_active());

  for (unsigned int j = 0; j < n2; ++j)
    for (unsigned int d = 0; d < nd; ++d)
      for (unsigned int k = 0; k < nq; ++k)
        matrix(k, j) += val2.shape_value_component(j, k, d);
}

int
main()
{
  const unsigned int fe_degree = 3;
  const double       factor    = 0.3;

  FE_DGQ<1> fe(fe_degree);
  QGauss<1> quad(fe_degree + 1);

  // compute projection matrix
  FullMatrix<double> P(fe.n_dofs_per_cell(), fe.n_dofs_per_cell());
  get_projection_matrix(fe, factor, P);
  P.print_formatted(std::cout);
  std::cout << std::endl;

  // compute composite interpolation matrix
  FullMatrix<double> I(quad.size(), fe.n_dofs_per_cell());
  get_interpolation_matrix(fe, quad, I);
  I.print_formatted(std::cout);
  std::cout << std::endl;

  FullMatrix<double> IxP(quad.size(), fe.n_dofs_per_cell());
  I.mmult(IxP, P);
  IxP.print_formatted(std::cout);
  std::cout << std::endl;

  // compute direct interpolation matrix
  get_interpolation_matrix(fe, scale(quad, factor), I);
  I.print_formatted(std::cout);
  std::cout << std::endl;

  // observaition: both interpolation matrices are the same
}
