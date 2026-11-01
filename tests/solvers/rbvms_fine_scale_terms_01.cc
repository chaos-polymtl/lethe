// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the RBVMS fine-scale terms of Bazilevs et al. 2007
 * (eq. 72) implemented by rbvms_fine_scale_terms (cross-stress +
 * Reynolds-stress gradient-residual contribution) and its consistent
 * linearization rbvms_fine_scale_terms_linearization. Both helpers are written
 * in terms of the fine-scale velocity u' and its increment δu'.
 *
 * Three checks are performed:
 *   1. With the quasi-static fine-scale velocity u' = -τ_M r_M, the residual
 *      contribution T[i][j] = -u_i u'_j - u'_i u'_j must equal the closed form
 *      τ_M u_i r_M[j] - τ_M^2 r_M[i] r_M[j] written in terms of the residual.
 *   2. With δu' = -τ_M δr_M, the linearization must equal the closed form
 *      τ_M (δu_i r_M[j] + u_i δr_M[j]) - τ_M^2 (δr_M[i] r_M[j] + r_M[i]
 *      δr_M[j]).
 *   3. The linearization is compared to a central finite difference of the
 *      residual contribution taken along (δu, δu'). Because the residual
 *      contribution is quadratic in its arguments, the central difference is
 *      exact up to round-off, so this is a strong consistency check of the
 *      Jacobian.
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
    const double tau_ref      = 0.25;
    const Number tau          = tau_ref;

    const Tensor<1, dim, Number> u   = make_vector<Number>(u0_arr);
    const Tensor<1, dim, Number> rm  = make_vector<Number>(rm0_arr);
    const Tensor<1, dim, Number> du  = make_vector<Number>(du_arr);
    const Tensor<1, dim, Number> drm = make_vector<Number>(drm_arr);

    // Quasi-static fine-scale velocity and its increment.
    const Tensor<1, dim, Number> up  = -tau * rm;
    const Tensor<1, dim, Number> dup = -tau * drm;

    // Tolerance for the comparison with the closed forms, which only differ by
    // the grouping of the floating-point operations.
    const double closed_form_tol = 1e-14;
    // Tolerance for the finite-difference check, which is exact up to the
    // round-off amplified by the finite-difference step.
    const double finite_difference_tol = 1e-9;

    // --- Check 1: residual contribution vs closed form in r_M. ---
    const Tensor<2, dim, Number> T = rbvms_fine_scale_terms<dim>(u, up);

    bool residual_ok = true;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        {
          const double ref = tau_ref * u0_arr[i] * rm0_arr[j] -
                             tau_ref * tau_ref * rm0_arr[i] * rm0_arr[j];
          if (std::abs(lane0(T[i][j]) - ref) > closed_form_tol)
            residual_ok = false;
        }
    deallog << label << " residual terms match the closed form in r_M      : "
            << (residual_ok ? "OK" : "FAIL") << std::endl;

    // --- Check 2: linearization vs closed form in (δu, δr_M). ---
    const Tensor<2, dim, Number> J =
      rbvms_fine_scale_terms_linearization<dim>(u, du, up, dup);

    bool linearization_ok = true;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        {
          const double ref =
            tau_ref * (du_arr[i] * rm0_arr[j] + u0_arr[i] * drm_arr[j]) -
            tau_ref * tau_ref *
              (drm_arr[i] * rm0_arr[j] + rm0_arr[i] * drm_arr[j]);
          if (std::abs(lane0(J[i][j]) - ref) > closed_form_tol)
            linearization_ok = false;
        }
    deallog << label << " linearization matches the closed form in r_M     : "
            << (linearization_ok ? "OK" : "FAIL") << std::endl;

    // --- Check 3: linearization vs central finite difference. ---
    const Number           h = 1e-2;
    Tensor<1, dim, Number> u_plus, u_minus, up_plus, up_minus;
    for (int i = 0; i < dim; ++i)
      {
        u_plus[i]   = u[i] + h * du[i];
        u_minus[i]  = u[i] - h * du[i];
        up_plus[i]  = up[i] + h * dup[i];
        up_minus[i] = up[i] - h * dup[i];
      }
    const Tensor<2, dim, Number> R_plus =
      rbvms_fine_scale_terms<dim>(u_plus, up_plus);
    const Tensor<2, dim, Number> R_minus =
      rbvms_fine_scale_terms<dim>(u_minus, up_minus);

    bool jacobian_ok = true;
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
        {
          const double fd =
            lane0((R_plus[i][j] - R_minus[i][j]) / (Number(2.0) * h));
          if (std::abs(lane0(J[i][j]) - fd) > finite_difference_tol)
            jacobian_ok = false;
        }
    deallog << label << " linearization matches finite difference          : "
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
