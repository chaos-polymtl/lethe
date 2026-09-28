// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the Gaussian and top-hat filter kernels used by the
 * Anderson-Jackson filter.
 *
 * For each kernel, in 2D and 3D, the test checks the value at the center of
 * the kernel, the compact support (the kernel is positive just inside the
 * support radius and zero on and beyond it) and the unit mass of the kernel,
 * integrated with a composite Gauss rule in the radial direction. For the
 * truncated Gaussian, the retained mass fraction given by the closed form is
 * also compared with its numerical integral. Finally, the test checks that
 * invalid kernels throw.
 */

// Deal.II includes
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>

// Lethe
#include <core/filter_kernels.h>

#include <cmath>
#include <iomanip>
#include <numbers>
#include <sstream>
#include <string>

// Tests (with common definitions)
#include <../tests/tests.h>

using namespace dealii;

namespace
{
  /**
   * @brief Format a number in scientific notation.
   *
   * @param[in] value Number to format.
   *
   * @param[in] precision Number of digits after the decimal point.
   *
   * @return Formatted number.
   */
  std::string
  format(const double value, const unsigned int precision = 6)
  {
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(precision) << value;
    return stream.str();
  }

  /**
   * @brief Format a boolean as true or false.
   *
   * @param[in] value Boolean to format.
   *
   * @return Formatted boolean.
   */
  std::string
  format(const bool value)
  {
    return value ? "true" : "false";
  }

  /**
   * @brief Integrate a radial function over a ball centered at the origin
   * with a composite Gauss rule in the radial direction.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @tparam FunctionType Type of the radial function.
   *
   * @param[in] radial_function Function of the distance to the origin.
   *
   * @param[in] radius Radius of the ball.
   *
   * @return Integral of the function over the ball.
   */
  template <int dim, typename FunctionType>
  double
  integrate_over_ball(const FunctionType &radial_function, const double radius)
  {
    const QIterated<1> quadrature(QGauss<1>(10), 50);

    double integral = 0.;
    for (unsigned int q = 0; q < quadrature.size(); ++q)
      {
        const double r = radius * quadrature.point(q)[0];
        // Measure of the sphere of radius r
        const double sphere_measure = (dim == 2) ?
                                        2. * std::numbers::pi * r :
                                        4. * std::numbers::pi * r * r;
        integral +=
          radial_function(r) * sphere_measure * quadrature.weight(q) * radius;
      }
    return integral;
  }

  /**
   * @brief Check the value at the center, the compact support and the unit
   * mass of a kernel.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @tparam KernelType Type of the kernel.
   *
   * @param[in] kernel Kernel to check.
   */
  template <int dim, typename KernelType>
  void
  check_kernel(const KernelType &kernel)
  {
    const double radius = kernel.support_radius();

    deallog << "  support radius                  : " << format(radius)
            << std::endl;
    deallog << "  squared support radius matches  : "
            << format(std::abs(kernel.support_radius_squared() -
                               radius * radius) <= 1e-14 * radius * radius)
            << std::endl;
    deallog << "  value at the center             : "
            << format(kernel.value_from_squared_distance(0.)) << std::endl;
    deallog << "  positive just inside support    : "
            << format(kernel.value_from_squared_distance((1. - 1e-10) * radius *
                                                         radius) > 0.)
            << std::endl;
    deallog << "  zero on the support radius      : "
            << format(kernel.value_from_squared_distance(radius * radius) == 0.)
            << std::endl;
    deallog << "  zero beyond the support radius  : "
            << format(
                 kernel.value_from_squared_distance(4. * radius * radius) == 0.)
            << std::endl;

    const double mass = integrate_over_ball<dim>(
      [&kernel](const double r) {
        return kernel.value_from_squared_distance(r * r);
      },
      radius);
    deallog << "  unit mass (|mass - 1| < 1e-12)  : "
            << format(std::abs(mass - 1.) < 1e-12) << std::endl;
  }

  /**
   * @brief Test the truncated Gaussian kernel.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] standard_deviation Standard deviation of the Gaussian.
   *
   * @param[in] cutoff Truncation radius in number of standard deviations.
   */
  template <int dim>
  void
  test_gaussian(const double standard_deviation, const double cutoff)
  {
    deallog << "Gaussian kernel, dim = " << dim
            << ", standard deviation = " << format(standard_deviation)
            << ", cutoff = " << format(cutoff) << std::endl;

    const GaussianFilterKernel<dim> kernel(standard_deviation, cutoff);
    check_kernel<dim>(kernel);

    deallog << "  retained mass fraction          : "
            << format(kernel.retained_mass_fraction()) << std::endl;

    // Mass of the untruncated Gaussian contained in the support, integrated
    // numerically.
    const double variance      = standard_deviation * standard_deviation;
    const double retained_mass = integrate_over_ball<dim>(
      [variance](const double r) {
        return std::exp(-0.5 * r * r / variance) /
               std::pow(2. * std::numbers::pi * variance, 0.5 * dim);
      },
      kernel.support_radius());
    deallog << "  retained mass fraction matches  : "
            << format(std::abs(retained_mass -
                               kernel.retained_mass_fraction()) < 1e-12)
            << std::endl;
  }

  /**
   * @brief Test the top-hat kernel.
   *
   * @tparam dim Number of spatial dimensions.
   *
   * @param[in] radius Radius of the kernel.
   */
  template <int dim>
  void
  test_top_hat(const double radius)
  {
    deallog << "Top-hat kernel, dim = " << dim
            << ", radius = " << format(radius) << std::endl;

    const TopHatFilterKernel<dim> kernel(radius);
    check_kernel<dim>(kernel);
  }

  /**
   * @brief Check that kernels with invalid parameters throw.
   *
   * @tparam dim Number of spatial dimensions.
   */
  template <int dim>
  void
  test_invalid_kernels()
  {
    deallog << "Invalid kernels, dim = " << dim << std::endl;

    bool zero_standard_deviation_throws = false;
    try
      {
        const GaussianFilterKernel<dim> kernel(0., 3.);
      }
    // The message of the exception is not printed since it contains the file
    // name and the line number, which are not portable across builds.
    catch (const ExceptionBase &)
      {
        zero_standard_deviation_throws = true;
      }
    deallog << "  zero standard deviation throws  : "
            << format(zero_standard_deviation_throws) << std::endl;

    bool zero_cutoff_throws = false;
    try
      {
        const GaussianFilterKernel<dim> kernel(1., 0.);
      }
    catch (const ExceptionBase &)
      {
        zero_cutoff_throws = true;
      }
    deallog << "  zero cutoff throws              : "
            << format(zero_cutoff_throws) << std::endl;

    bool negative_radius_throws = false;
    try
      {
        const TopHatFilterKernel<dim> kernel(-1.);
      }
    catch (const ExceptionBase &)
      {
        negative_radius_throws = true;
      }
    deallog << "  negative top-hat radius throws  : "
            << format(negative_radius_throws) << std::endl;
  }
} // namespace

int
main()
{
  try
    {
      initlog();

      test_gaussian<2>(0.5, 3.);
      test_gaussian<3>(0.5, 3.);
      test_gaussian<2>(2., 1.5);
      test_gaussian<3>(2., 1.5);
      test_top_hat<2>(0.7);
      test_top_hat<3>(0.7);
      test_invalid_kernels<2>();
      test_invalid_kernels<3>();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl
                << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}
