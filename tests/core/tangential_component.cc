// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Test the tangential_component function in 2D and 3D, for real and
 * complex vectors, by comparing it with the exact tangential component
 * v - (v.n)n of the vector.
 */

// Deal.II
#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

// Lethe
#include <core/tensors_and_points_dimension_manipulation.h>

// Tests
#include <../tests/tests.h>

#include <cmath>
#include <complex>
#include <string>

/**
 * @brief Compute the tangential component of a vector and check that it is
 * equal to the expected one.
 *
 * @param[in] name Name of the case printed in the output.
 * @param[in] tensor Vector whose tangential component is computed.
 * @param[in] normal Unit normal vector of the surface.
 * @param[in] expected Exact tangential component of the vector.
 */
template <int dim, typename Number>
void
check_case(const std::string            &name,
           const Tensor<1, dim, Number> &tensor,
           const Tensor<1, dim>         &normal,
           const Tensor<1, dim, Number> &expected)
{
  const double tolerance = 1e-12;

  const Tensor<1, dim, Number> result = tangential_component(tensor, normal);

  AssertThrow((result - expected).norm() < tolerance,
              ExcMessage("The tangential component of the case " + name +
                         " is not equal to the expected one."));
  deallog << name << ": OK" << std::endl;
}

void
test()
{
  const std::complex<double> imag(0., 1.);

  // 2D, real vector
  {
    const Tensor<1, 2> v({3., 4.});
    check_case("2D real, n = (1, 0)",
               v,
               Tensor<1, 2>({1., 0.}),
               Tensor<1, 2>({0., 4.}));
    check_case("2D real, n = (1, 1)/sqrt(2)",
               v,
               Tensor<1, 2>({1. / std::sqrt(2.), 1. / std::sqrt(2.)}),
               Tensor<1, 2>({-0.5, 0.5}));
  }

  // 2D, complex vector
  {
    const Tensor<1, 2, std::complex<double>> v({1. + 2. * imag, 3. - imag});
    check_case("2D complex, n = (0, 1)",
               v,
               Tensor<1, 2>({0., 1.}),
               Tensor<1, 2, std::complex<double>>({1. + 2. * imag, 0.}));
  }

  // 3D, real vector
  {
    const Tensor<1, 3> v({1., 2., 3.});
    check_case("3D real, n = (0, 0, 1)",
               v,
               Tensor<1, 3>({0., 0., 1.}),
               Tensor<1, 3>({1., 2., 0.}));
    check_case("3D real, n = (1, 1, 1)/sqrt(3)",
               v,
               Tensor<1, 3>(
                 {1. / std::sqrt(3.), 1. / std::sqrt(3.), 1. / std::sqrt(3.)}),
               Tensor<1, 3>({-1., 0., 1.}));
  }

  // 3D, complex vector
  {
    const Tensor<1, 3, std::complex<double>> v(
      {1. + imag, 2. - imag, 3. * imag});
    check_case("3D complex, n = (0, 0, 1)",
               v,
               Tensor<1, 3>({0., 0., 1.}),
               Tensor<1, 3, std::complex<double>>({1. + imag, 2. - imag, 0.}));
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
                << "Aborting!" << std::endl;
      return 1;
    }

  return 0;
}
