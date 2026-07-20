// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Matrix-based counterpart of rbvms_metric_tensor. It verifies that the
 * RBVMS metric tensor G (eq. 66) and vector g (eq. 69) are computed correctly
 * when the inverse Jacobian is obtained from FEValues rather than from the
 * matrix-free FEEvaluation.
 *
 * The subtlety being tested is the storage convention:
 * FEValues::inverse_jacobian returns [i][j] = dxi_i/dx_j, which is the
 * TRANSPOSE of the matrix-free FEEvaluation::inverse_jacobian ([i][j] =
 * dxi_j/dx_i) that compute_metric_tensor expects. The matrix-based RBVMS
 * assembler therefore transposes the FEValues inverse Jacobian before calling
 * compute_metric_tensor; this test replicates and checks that handling on the
 * same stretch-and-shear cell used by rbvms_metric_tensor (x = A*xi with A =
 * [[2,1,0],[0,2,1],[0,0,2]]), so the analytic G, g and the tau-scalars are
 * identical and checked to 1e-12.
 */

// Deal.II
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

// Lethe (provides compute_metric_tensor)
#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

void
test()
{
  const int dim = 3;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  GridTools::transform(
    [](const Point<dim> &p) {
      return Point<dim>(2. * p[0] + 1. * p[1],
                        2. * p[1] + 1. * p[2],
                        2. * p[2]);
    },
    tria);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  MappingQ1<dim> mapping;
  QGauss<dim>    quadrature(2);

  FEValues<dim> fe_values(mapping, fe, quadrature, update_inverse_jacobians);
  fe_values.reinit(dof_handler.begin_active());

  // FEValues inverse Jacobian: [i][j] = dxi_i/dx_j. Transpose it into the
  // matrix-free convention ([i][j] = dxi_j/dx_i) expected by
  // compute_metric_tensor, exactly as the RBVMS assembler does.
  const auto            &fe_inverse_jacobian = fe_values.inverse_jacobian(0);
  Tensor<2, dim, double> inverse_jacobian;
  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < dim; ++j)
      inverse_jacobian[i][j] = fe_inverse_jacobian[j][i];

  Tensor<2, dim, double> G;
  Tensor<1, dim, double> g;
  compute_metric_tensor(inverse_jacobian, G, g);

  // Analytic reference values (identical to rbvms_metric_tensor).
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

  deallog << "metric tensor G (eq. 66) matches analytic : "
          << ((G - G_exact).norm() < tol ? "OK" : "FAIL") << std::endl;
  deallog << "metric vector g (eq. 69) matches analytic  : "
          << ((g - g_exact).norm() < tol ? "OK" : "FAIL") << std::endl;

  Tensor<1, dim, double> u;
  u[0] = 1.;
  u[1] = 2.;
  u[2] = 3.;

  double uGu = 0., GG = 0., gg = 0.;
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
