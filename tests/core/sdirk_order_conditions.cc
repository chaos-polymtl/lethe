// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Verify that the Butcher tableaux of the SDIRK methods satisfy the
 * Runge-Kutta order conditions up to the order of each method and that the
 * stage times c are the row sums of A. The test also reports whether each
 * method is stiffly accurate (b is the last row of A) and the value of its
 * stability function at infinity, R(inf) = 1 - b^T A^{-1} 1, which vanishes
 * for L-stable stiffly accurate methods.
 */

// Lethe
#include <core/parameters.h>
#include <core/sdirk_stage_data.h>

// deal.II
#include <deal.II/lac/full_matrix.h>

// Tests (with common definitions)
#include <../tests/tests.h>

// Std
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>

namespace
{
  /// Tolerance on the order conditions. The coefficients of the tableaux are
  /// expected to be accurate to machine precision.
  constexpr double order_condition_tolerance = 1e-13;

  /**
   * @brief Dot product of two vectors of the same size.
   *
   * @param[in] x First vector.
   * @param[in] y Second vector.
   *
   * @return The dot product x^T y.
   */
  double
  dot(const std::vector<double> &x, const std::vector<double> &y)
  {
    double result = 0.;
    for (unsigned int i = 0; i < x.size(); ++i)
      result += x[i] * y[i];
    return result;
  }

  /**
   * @brief Element-wise product of two vectors of the same size.
   *
   * @param[in] x First vector.
   * @param[in] y Second vector.
   *
   * @return The vector whose entries are x_i * y_i.
   */
  std::vector<double>
  elementwise_product(const std::vector<double> &x,
                      const std::vector<double> &y)
  {
    std::vector<double> result(x.size());
    for (unsigned int i = 0; i < x.size(); ++i)
      result[i] = x[i] * y[i];
    return result;
  }

  /**
   * @brief Matrix-vector product.
   *
   * @param[in] A Matrix.
   * @param[in] x Vector.
   *
   * @return The vector A x.
   */
  std::vector<double>
  matrix_vector_product(const FullMatrix<double>  &A,
                        const std::vector<double> &x)
  {
    std::vector<double> result(A.m(), 0.);
    for (unsigned int i = 0; i < A.m(); ++i)
      for (unsigned int j = 0; j < A.n(); ++j)
        result[i] += A(i, j) * x[j];
    return result;
  }

  /**
   * @brief Residuals of the Runge-Kutta order conditions of order exactly p.
   *
   * @param[in] table Butcher tableau of the method.
   * @param[in] p Order of the conditions (1 to 4).
   *
   * @return One residual per rooted tree of order p.
   */
  std::vector<double>
  order_condition_residuals(const SDIRKTable &table, const unsigned int p)
  {
    const FullMatrix<double>  &A = table.A;
    const std::vector<double> &b = table.b;
    const std::vector<double> &c = table.c;

    const std::vector<double> ones(b.size(), 1.);
    const std::vector<double> c2 = elementwise_product(c, c);
    const std::vector<double> Ac = matrix_vector_product(A, c);

    switch (p)
      {
        case 1:
          return {dot(b, ones) - 1.};
        case 2:
          return {dot(b, c) - 1. / 2.};
        case 3:
          return {dot(b, c2) - 1. / 3., dot(b, Ac) - 1. / 6.};
        case 4:
          return {dot(b, elementwise_product(c2, c)) - 1. / 4.,
                  dot(b, elementwise_product(c, Ac)) - 1. / 8.,
                  dot(b, matrix_vector_product(A, c2)) - 1. / 12.,
                  dot(b, matrix_vector_product(A, Ac)) - 1. / 24.};
        default:
          AssertThrow(false,
                      ExcMessage("The order conditions are only implemented "
                                 "up to order 4."));
          return {};
      }
  }

  /**
   * @brief Value of the stability function at infinity,
   * R(inf) = 1 - b^T A^{-1} 1, obtained by forward substitution since A is
   * lower triangular.
   *
   * @param[in] table Butcher tableau of the method.
   *
   * @return R(inf).
   */
  double
  stability_function_at_infinity(const SDIRKTable &table)
  {
    const unsigned int  n_stages = table.b.size();
    std::vector<double> x(n_stages);
    for (unsigned int i = 0; i < n_stages; ++i)
      {
        double rhs = 1.;
        for (unsigned int j = 0; j < i; ++j)
          rhs -= table.A(i, j) * x[j];
        x[i] = rhs / table.A(i, i);
      }
    return 1. - dot(table.b, x);
  }
} // namespace

void
test()
{
  using Method = Parameters::SimulationControl::TimeSteppingMethod;

  deallog << std::setprecision(6) << std::fixed;

  struct MethodToTest
  {
    std::string  name;
    Method       method;
    unsigned int order;
  };

  // The naming convention of the SDIRK methods is sdirkOrderStages
  const std::vector<MethodToTest> methods = {{"sdirk22", Method::sdirk22, 2},
                                             {"sdirk33", Method::sdirk33, 3},
                                             {"sdirk45", Method::sdirk45, 4}};

  for (const auto &[name, method, order] : methods)
    {
      const SDIRKTable   table    = sdirk_table(method);
      const unsigned int n_stages = table.b.size();

      // The stage times are the row sums of A
      for (unsigned int i = 0; i < n_stages; ++i)
        {
          double row_sum = 0.;
          for (unsigned int j = 0; j < n_stages; ++j)
            row_sum += table.A(i, j);
          AssertThrow(std::abs(row_sum - table.c[i]) <
                        order_condition_tolerance,
                      ExcMessage("The c coefficients of " + name +
                                 " are not the row sums of A."));
        }

      // Order conditions up to the order of the method
      for (unsigned int p = 1; p <= order; ++p)
        for (const double residual : order_condition_residuals(table, p))
          AssertThrow(std::abs(residual) < order_condition_tolerance,
                      ExcMessage("An order condition of order " +
                                 std::to_string(p) + " is not satisfied by " +
                                 name + "."));

      // A method is stiffly accurate if b is the last row of A
      bool stiffly_accurate = true;
      for (unsigned int j = 0; j < n_stages; ++j)
        if (std::abs(table.A(n_stages - 1, j) - table.b[j]) >
            order_condition_tolerance)
          stiffly_accurate = false;

      double r_infinity = stability_function_at_infinity(table);
      if (std::abs(r_infinity) < order_condition_tolerance)
        r_infinity = 0.;

      deallog << name << ": " << n_stages
              << " stages, order conditions satisfied up to order " << order
              << std::endl;
      deallog << "  stiffly accurate: " << (stiffly_accurate ? "true" : "false")
              << std::endl;
      deallog << "  R(infinity): " << r_infinity << std::endl;
    }
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
