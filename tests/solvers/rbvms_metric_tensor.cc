// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the RBVMS element metric tensor G (Bazilevs et al.
 * 2007, eq. 66) and metric vector g (eq. 69) computed by compute_metric_tensor
 * from deal.II's FEEvaluation::inverse_jacobian(q), on a cell that is BOTH
 * stretched AND sheared. Because G is symmetric, a wrong transpose convention
 * on inverse_jacobian only shows up on a sheared cell, so the sheared cell is
 * the discriminating test. The physical cell is obtained from the unit cube by
 * the affine map x = A*xi with the upper-triangular matrix
 *
 *        [ 2 1 0 ]
 *    A = [ 0 2 1 ]
 *        [ 0 0 2 ]
 *
 * so the mapping Jacobian is J = A (constant, affine cell) and inverse_jacobian
 * returns J^{-T}. The analytic values of G, g and the tau-scalars (u.G.u, G:G,
 * g.g) for u = (1,2,3) are hard-coded below and checked to 1e-12.
 */

// Deal.II
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

// Lethe (provides compute_metric_tensor)
#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  const int dim = 3;

  // Build the unit cube and apply the stretch-and-shear affine map x = A*xi.
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  GridTools::transform(
    [](const Point<dim> &p) {
      return Point<dim>(2. * p[0] + 1. * p[1],
                        2. * p[1] + 1. * p[2],
                        2. * p[2]);
    },
    tria);

  // Minimal matrix-free setup: a scalar FE is sufficient since inverse_jacobian
  // is a purely geometric quantity.
  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<double> constraints;
  constraints.close();

  MappingQ1<dim> mapping;
  QGauss<1>      quadrature(2);

  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.mapping_update_flags =
    update_gradients | update_JxW_values | update_quadrature_points;

  MatrixFree<dim, double> matrix_free;
  matrix_free.reinit(
    mapping, dof_handler, constraints, quadrature, additional_data);

  FEEvaluation<dim, -1, 0, 1, double> fe_eval(matrix_free);
  fe_eval.reinit(0);

  // Read deal.II's inverse Jacobian (J^{-T}) at the first quadrature point and
  // extract lane 0 (single-cell batch).
  const auto             ij_vec = fe_eval.inverse_jacobian(0);
  Tensor<2, dim, double> inv_j;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      inv_j[i][j] = ij_vec[i][j][0];

  // Compute the metric tensor and vector using the helper under test.
  Tensor<2, dim, double> G;
  Tensor<1, dim, double> g;
  compute_metric_tensor(inv_j, G, g);

  // Analytic reference values (see file header), exact binary fractions.
  Tensor<2, dim, double> G_exact;
  G_exact[0][0] = 0.25;
  G_exact[0][1] = -0.125;
  G_exact[0][2] = 0.0625;
  G_exact[1][0] = -0.125;
  G_exact[1][1] = 0.3125;
  G_exact[1][2] = -0.15625;
  G_exact[2][0] = 0.0625;
  G_exact[2][1] = -0.15625;
  G_exact[2][2] = 0.328125;

  Tensor<1, dim, double> g_exact;
  g_exact[0] = 0.5;
  g_exact[1] = 0.25;
  g_exact[2] = 0.375;

  const double tol = 1e-12;

  // Check G (eq. 66).
  double g_tensor_error = (G - G_exact).norm();
  deallog << "metric tensor G (eq. 66) matches analytic : "
          << (g_tensor_error < tol ? "OK" : "FAIL") << std::endl;

  // Check symmetry of G explicitly (eq. 66 states G is symmetric).
  double g_symmetry_error = (G - transpose(G)).norm();
  deallog << "metric tensor G is symmetric               : "
          << (g_symmetry_error < tol ? "OK" : "FAIL") << std::endl;

  // Check g (eq. 69).
  double g_vector_error = (g - g_exact).norm();
  deallog << "metric vector g (eq. 69) matches analytic  : "
          << (g_vector_error < tol ? "OK" : "FAIL") << std::endl;

  // Check the tau scalars (eqs. 67, 68, 70) for u = (1, 2, 3).
  Tensor<1, dim, double> u;
  u[0] = 1.;
  u[1] = 2.;
  u[2] = 3.;

  double uGu = 0.; // eq. 68
  double GG  = 0.; // eq. 67
  double gg  = 0.; // eq. 70
  for (int i = 0; i < dim; ++i)
    {
      gg += g[i] * g[i];
      for (int j = 0; j < dim; ++j)
        {
          uGu += u[i] * G[i][j] * u[j];
          GG += G[i][j] * G[i][j];
        }
    }

  deallog << "u.G.u (eq. 68) matches analytic            : "
          << (std::abs(uGu - 2.453125) < tol ? "OK" : "FAIL") << std::endl;
  deallog << "G:G   (eq. 67) matches analytic            : "
          << (std::abs(GG - 0.355712890625) < tol ? "OK" : "FAIL") << std::endl;
  deallog << "g.g   (eq. 70) matches analytic            : "
          << (std::abs(gg - 0.453125) < tol ? "OK" : "FAIL") << std::endl;

  // Transpose sentinel: feeding the transposed inverse Jacobian must yield a
  // DIFFERENT G on this sheared cell. This guards against the deal.II storage
  // convention flipping in a future version (see compute_metric_tensor docs).
  Tensor<2, dim, double> G_transposed;
  Tensor<1, dim, double> g_transposed;
  compute_metric_tensor(transpose(inv_j), G_transposed, g_transposed);
  double sentinel_diff = (G - G_transposed).norm();
  deallog << "transpose sentinel (G is shear-sensitive)  : "
          << (sentinel_diff > 1e-6 ? "OK" : "FAIL") << std::endl;
}

int
main(int argc, char **argv)
{
  try
    {
      initlog();
      Utilities::MPI::MPI_InitFinalize mpi_initialization(
        argc, argv, numbers::invalid_unsigned_int);
      test();
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
