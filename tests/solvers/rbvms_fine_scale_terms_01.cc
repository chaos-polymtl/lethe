// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the RBVMS fine-scale terms of Bazilevs et al. 2007
 * (eq. 72) implemented by rbvms_fine_scale_terms (cross-stress +
 * Reynolds-stress gradient-residual contribution) and its consistent frozen-τ
 * linearization rbvms_fine_scale_terms_linearization.
 *
 * Two checks are performed:
 *   1. The residual contribution T[i][j] = τ_M u_i r_M[j] - τ_M^2 r_M[i] r_M[j]
 *      is compared entry-by-entry to an independent in-test evaluation.
 *   2. The linearization is compared to a central finite difference of the
 *      residual contribution taken along (δu, δr_M) with τ_M held fixed.
 * Because the residual contribution is quadratic in its arguments, the central
 *      difference is exact up to round-off, so this is a strong consistency
 *      check of the Jacobian (the frozen-τ derivative of the residual).
 *
 * Both double and VectorizedArray<double> instantiations are exercised.
 */

// Deal.II
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

// Lethe
#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

namespace
{
  constexpr int dim = 3;

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

  // Fill a dim-vector from a double array (works for double and vectorized).
  template <typename Number>
  Tensor<1, dim, Number>
  make_vector(const double (&a)[dim])
  {
    Tensor<1, dim, Number> v;
    for (int i = 0; i < dim; ++i)
      v[i] = a[i];
    return v;
  }

  template <typename Number>
  void
  run(const std::string &label)
  {
    const double u0_arr[dim]  = {1.0, 2.0, 3.0};
    const double rm0_arr[dim] = {0.5, -1.0, 2.0};
    const double du_arr[dim]  = {0.1, -0.2, 0.3};
    const double drm_arr[dim] = {-0.5, 0.4, 0.7};
    const Number tau          = 0.25;

    const Tensor<1, dim, Number> u   = make_vector<Number>(u0_arr);
    const Tensor<1, dim, Number> rm  = make_vector<Number>(rm0_arr);
    const Tensor<1, dim, Number> du  = make_vector<Number>(du_arr);
    const Tensor<1, dim, Number> drm = make_vector<Number>(drm_arr);

    const double tol = 1e-9;

    // --- Check 1: residual contribution vs independent reference. ---
    const Tensor<2, dim, Number> T = rbvms_fine_scale_terms<dim>(u, rm, tau);

    const double tau_ref     = 0.25;
    bool         residual_ok = true;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        {
          const double ref = tau_ref * u0_arr[i] * rm0_arr[j] -
                             tau_ref * tau_ref * rm0_arr[i] * rm0_arr[j];
          if (std::abs(lane0(T[i][j]) - ref) > tol)
            residual_ok = false;
        }
    deallog << label << " residual terms (eq. 72) match reference : "
            << (residual_ok ? "OK" : "FAIL") << std::endl;

    // --- Check 2: linearization vs central finite difference. ---
    const Tensor<2, dim, Number> J =
      rbvms_fine_scale_terms_linearization<dim>(u, du, rm, drm, tau);

    const Number           h = 1e-2;
    Tensor<1, dim, Number> u_plus, u_minus, rm_plus, rm_minus;
    for (int i = 0; i < dim; ++i)
      {
        u_plus[i]   = u[i] + h * du[i];
        u_minus[i]  = u[i] - h * du[i];
        rm_plus[i]  = rm[i] + h * drm[i];
        rm_minus[i] = rm[i] - h * drm[i];
      }
    const Tensor<2, dim, Number> R_plus =
      rbvms_fine_scale_terms<dim>(u_plus, rm_plus, tau);
    const Tensor<2, dim, Number> R_minus =
      rbvms_fine_scale_terms<dim>(u_minus, rm_minus, tau);

    bool jacobian_ok = true;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        {
          const double fd =
            lane0((R_plus[i][j] - R_minus[i][j]) / (Number(2.0) * h));
          if (std::abs(lane0(J[i][j]) - fd) > tol)
            jacobian_ok = false;
        }
    deallog << label << " linearization matches finite difference  : "
            << (jacobian_ok ? "OK" : "FAIL") << std::endl;
  }
} // namespace

void
test()
{
  run<double>("[double]    ");
  run<VectorizedArray<double>>("[vectorized]");
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
