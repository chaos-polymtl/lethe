// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test verifies the discrete dynamic RBVMS fine-scale velocity
 *
 *   u'_{n+1} = -tau_t (r_M + m'),   tau_t = (beta_0 + 1/tau_s)^(-1),
 *   m' = sum_{j>=1} beta_j u'_{n+1-j},
 *
 * which is the BDF discretization of d(u')/dt + u'/tau_s = -r_M. The
 * stabilization parameter tau_t is computed with calculate_rbvms_dynamic_tau
 * and the BDF coefficients beta_j with calculate_bdf_coefficients.
 *
 * Three checks are performed:
 *   1. In the steady limit (alpha = 0), tau_t is equal to tau_s.
 *   2. For BDF1, u'_{n+1} = tau_t (-r_M + u'_n/dt), i.e. the backward Euler
 *      discretization of the fine-scale equation.
 *   3. For BDF1, BDF2 and BDF3 with a constant residual and a constant time
 *      step, the recursion converges to the steady fine-scale velocity
 *      -tau_s r_M, whatever the time step. This holds because the BDF
 *      coefficients sum to zero, and it is the property that makes the steady
 *      solution independent of the time step.
 *
 * Checks 1 and 2 are performed for both double and VectorizedArray<double>.
 */

// Deal.II
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/lac/vector.h>

// Lethe
#include <core/bdf.h>
#include <core/parameters.h>

#include <solvers/stabilization.h>

// Tests
#include <../tests/tests.h>

// C++
#include <array>
#include <string>
#include <vector>

namespace
{
  constexpr int dim = 2;

  using Method = Parameters::SimulationControl::TimeSteppingMethod;

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

  // Steady stabilization parameter and constant strong momentum residual used
  // by all the checks.
  constexpr double tau_static_value = 0.25;
  constexpr double r_m_values[dim]  = {0.5, -1.5};

  template <typename Number>
  void
  check_pointwise(const std::string &label)
  {
    // Tolerance for quantities that only differ by the grouping of the
    // floating-point operations.
    const double tolerance = 1e-14;

    const Number tau_static = tau_static_value;

    // --- Check 1: steady limit. ---
    const Number tau_steady = calculate_rbvms_dynamic_tau(tau_static, 0.);
    deallog << label << " steady limit tau_t = tau_s          : "
            << (std::abs(lane0(tau_steady) - tau_static_value) < tolerance ?
                  "OK" :
                  "FAIL")
            << std::endl;

    // --- Check 2: BDF1 is the backward Euler discretization. ---
    const double              dt = 0.1;
    const std::vector<double> time_steps(1, dt);
    const Vector<double>      bdf_coefs =
      calculate_bdf_coefficients(Method::bdf1, time_steps);

    Tensor<1, dim, Number> r_m;
    Tensor<1, dim, Number> previous_fine_scale_velocity;
    for (int d = 0; d < dim; ++d)
      {
        r_m[d]                          = r_m_values[d];
        previous_fine_scale_velocity[d] = 0.3 * (d + 1);
      }

    const Number tau = calculate_rbvms_dynamic_tau(tau_static, bdf_coefs[0]);
    const Tensor<1, dim, Number> memory =
      Number(bdf_coefs[1]) * previous_fine_scale_velocity;
    const Tensor<1, dim, Number> fine_scale_velocity = -tau * (r_m + memory);

    // Backward Euler: (u' - u'_n)/dt + u'/tau_s = -r_M
    const double tau_backward_euler = 1. / (1. / dt + 1. / tau_static_value);
    bool         backward_euler_ok  = true;
    for (int d = 0; d < dim; ++d)
      {
        const double reference =
          tau_backward_euler *
          (-r_m_values[d] + lane0(previous_fine_scale_velocity[d]) / dt);
        if (std::abs(lane0(fine_scale_velocity[d]) - reference) > tolerance)
          backward_euler_ok = false;
      }
    deallog << label << " BDF1 matches backward Euler         : "
            << (backward_euler_ok ? "OK" : "FAIL") << std::endl;
  }

  void
  check_steady_state(const Method method, const double dt)
  {
    // Tolerance on the converged fine-scale velocity.
    const double tolerance = 1e-12;
    // Enough steps for the slowest case (dt/tau_s = 0.04) to converge to
    // round-off.
    const unsigned int n_steps = 5000;

    const std::vector<double> time_steps(maximum_number_of_previous_solutions(),
                                         dt);
    const Vector<double>      bdf_coefs =
      calculate_bdf_coefficients(method, time_steps);
    const unsigned int n_previous = number_of_previous_solutions(method);

    const double tau =
      calculate_rbvms_dynamic_tau(tau_static_value, bdf_coefs[0]);

    Tensor<1, dim> r_m;
    for (int d = 0; d < dim; ++d)
      r_m[d] = r_m_values[d];

    // history[p] holds u'_{n-p}; the fine-scale velocity starts from zero.
    std::array<Tensor<1, dim>, 3> history;
    for (unsigned int step = 0; step < n_steps; ++step)
      {
        Tensor<1, dim> memory;
        for (unsigned int p = 0; p < n_previous; ++p)
          memory += bdf_coefs[p + 1] * history[p];

        const Tensor<1, dim> fine_scale_velocity = -tau * (r_m + memory);

        history[2] = history[1];
        history[1] = history[0];
        history[0] = fine_scale_velocity;
      }

    const Tensor<1, dim> steady_fine_scale_velocity = -tau_static_value * r_m;
    deallog << "bdf" << n_previous << " dt = " << dt
            << " converges to -tau_s r_M : "
            << ((history[0] - steady_fine_scale_velocity).norm() < tolerance ?
                  "OK" :
                  "FAIL")
            << std::endl;
  }
} // namespace

void
test()
{
  check_pointwise<double>("[double]    ");
  check_pointwise<VectorizedArray<double>>("[vectorized]");

  for (const Method method : {Method::bdf1, Method::bdf2, Method::bdf3})
    for (const double dt : {0.01, 0.1, 1.0})
      check_steady_state(method, dt);
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
