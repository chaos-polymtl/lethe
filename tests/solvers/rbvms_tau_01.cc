// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the RBVMS stabilization parameters tau_M (eq. 64)
 * and tau_C (eq. 65) of Bazilevs et al. 2007 as implemented by
 * calculate_rbvms_tau. The metric tensor G and vector g are the analytic values
 * of the stretch-and-shear cell used in rbvms_metric_tensor. The function
 * output is checked against an independent in-test evaluation of eqs. 64-65 and
 * against the defining identity tau_C * tau_M * (g.g) = 1. The steady limit
 * (4/dt^2 = 0) is exercised as well. Both double and VectorizedArray<double>
 * instantiations are tested so the SIMD path is covered.
 */

// Deal.II
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

// Lethe (provides calculate_rbvms_tau)
#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

namespace
{
  constexpr int dim = 3;

  // Analytic metric tensor / vector of the stretch-and-shear cell (see
  // rbvms_metric_tensor).
  template <typename Number>
  void
  fill_metric(Tensor<2, dim, Number> &G, Tensor<1, dim, Number> &g)
  {
    G[0][0] = 0.25;
    G[0][1] = -0.125;
    G[0][2] = 0.0625;
    G[1][0] = -0.125;
    G[1][1] = 0.3125;
    G[1][2] = -0.15625;
    G[2][0] = 0.0625;
    G[2][1] = -0.15625;
    G[2][2] = 0.328125;

    g[0] = 0.5;
    g[1] = 0.25;
    g[2] = 0.375;
  }

  // Lane-0 extraction that works for both double and VectorizedArray<double>.
  double
  lane0(const double &x)
  {
    return x;
  }
  double
  lane0(const VectorizedArray<double> &x)
  {
    return x[0];
  }

  template <typename Number>
  void
  run(const std::string &label,
      const Number       four_over_dt_squared,
      const double       four_over_dt_squared_ref)
  {
    Tensor<2, dim, Number> G;
    Tensor<1, dim, Number> g;
    fill_metric(G, g);

    Tensor<1, dim, Number> u;
    u[0] = 1.0;
    u[1] = 2.0;
    u[2] = 3.0;

    const Number nu  = 0.5;
    const double c_i = 12.0; // = 3*k^2 for k = 2

    Number tau_m, tau_c;
    calculate_rbvms_tau<dim, Number>(
      G, g, u, nu, four_over_dt_squared, c_i, tau_m, tau_c);

    // Independent reference evaluation of eqs. 64-65 (double arithmetic).
    double G_arr[dim][dim];
    double g_arr[dim];
    double u_arr[dim];
    for (int i = 0; i < dim; ++i)
      {
        g_arr[i] = lane0(g[i]);
        u_arr[i] = lane0(u[i]);
        for (int j = 0; j < dim; ++j)
          G_arr[i][j] = lane0(G[i][j]);
      }
    double uGu = 0., GG = 0., gg = 0.;
    for (int i = 0; i < dim; ++i)
      {
        gg += g_arr[i] * g_arr[i];
        for (int j = 0; j < dim; ++j)
          {
            uGu += u_arr[i] * G_arr[i][j] * u_arr[j];
            GG += G_arr[i][j] * G_arr[i][j];
          }
      }
    const double nu_ref    = 0.5;
    const double tau_m_ref = 1.0 / std::sqrt(four_over_dt_squared_ref + uGu +
                                             c_i * nu_ref * nu_ref * GG);
    const double tau_c_ref = 1.0 / (tau_m_ref * gg);

    const double tol = 1e-12;

    deallog << label << " tau_M (eq. 64) matches reference : "
            << (std::abs(lane0(tau_m) - tau_m_ref) < tol ? "OK" : "FAIL")
            << std::endl;
    deallog << label << " tau_C (eq. 65) matches reference : "
            << (std::abs(lane0(tau_c) - tau_c_ref) < tol ? "OK" : "FAIL")
            << std::endl;
    // Defining identity tau_C * tau_M * (g.g) = 1 (eq. 65).
    deallog << label << " tau_C*tau_M*(g.g) = 1 identity   : "
            << (std::abs(lane0(tau_c) * lane0(tau_m) * gg - 1.0) < tol ? "OK" :
                                                                         "FAIL")
            << std::endl;
  }
} // namespace

void
test()
{
  // Transient: dt = 2 -> 4/dt^2 = 1.
  run<double>("[double,transient]     ", 1.0, 1.0);
  run<VectorizedArray<double>>("[vectorized,transient] ",
                               VectorizedArray<double>(1.0),
                               1.0);
  // Steady: 4/dt^2 = 0.
  run<double>("[double,steady]        ", 0.0, 0.0);
}

int
main()
{
  try
    {
      initlog();
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
