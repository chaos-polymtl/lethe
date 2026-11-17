// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief This test checks the stencil of the sharp-interface immersed boundary
 * method. The coefficients of the stencil must extrapolate exactly, to the
 * surface of the particle, the polynomials whose degree does not exceed the
 * degree of the stencil, and the derivative coefficients must give exactly the
 * derivatives of these polynomials on the surface. The interpolation points
 * must lie on the normal to the surface that passes through the degree of
 * freedom, and the derivative coefficients must give the derivative along this
 * normal of a quadratic function.
 */

// Lethe
#include <core/ib_particle.h>
#include <core/ib_stencil.h>

// Tests (with common definitions)
#include <../tests/tests.h>

#include <cmath>
#include <vector>

void
test()
{
  IBStencil<2> stencil;

  // Tolerance of the extrapolation of the monomials, relative to the sum of
  // the absolute values of the terms of the extrapolation. The coefficients
  // alternate in sign and grow with the degree and the length ratio, so the
  // round-off error scales with this sum rather than with the exact value.
  constexpr double extrapolation_tolerance = 1e-10;

  for (const double length_ratio : {1., 4., 8.})
    for (unsigned int degree = 1; degree <= 4; ++degree)
      {
        const std::vector<double> coefficients =
          stencil.coefficients(degree, length_ratio);

        // Coordinates of the points of the stencil on the reference line, in
        // the order of the coefficients: the degree of freedom (at 1) first,
        // then the interpolation points by increasing distance to the degree
        // of freedom. The surface of the particle is at 1 + length_ratio.
        std::vector<double> reference_points(degree + 1);
        for (unsigned int k = 0; k <= degree; ++k)
          reference_points[k] =
            0.5 * (1. - std::cos(numbers::PI * (degree - k) / degree));

        const std::vector<double> derivative_coefficients =
          stencil.derivative_coefficients(degree, length_ratio);

        // The monomial s^n and its derivative n s^(n-1) at the surface
        double maximal_relative_error            = 0.;
        double maximal_derivative_relative_error = 0.;
        for (unsigned int n = 0; n <= degree; ++n)
          {
            double extrapolated_value                = 0.;
            double sum_of_absolute_values            = 0.;
            double extrapolated_derivative           = 0.;
            double sum_of_absolute_derivative_values = 0.;
            for (unsigned int k = 0; k <= degree; ++k)
              {
                const double monomial = std::pow(reference_points[k], n);
                const double term     = coefficients[k] * monomial;
                extrapolated_value += term;
                sum_of_absolute_values += std::abs(term);

                const double derivative_term =
                  derivative_coefficients[k] * monomial;
                extrapolated_derivative += derivative_term;
                sum_of_absolute_derivative_values += std::abs(derivative_term);
              }
            const double exact_value = std::pow(1. + length_ratio, n);
            maximal_relative_error =
              std::max(maximal_relative_error,
                       std::abs(extrapolated_value - exact_value) /
                         sum_of_absolute_values);

            // The derivative of a constant is zero, so the error is only
            // relative to the sum of the absolute values of the terms
            const double exact_derivative =
              (n == 0) ? 0. : n * std::pow(1. + length_ratio, n - 1);
            maximal_derivative_relative_error =
              std::max(maximal_derivative_relative_error,
                       std::abs(extrapolated_derivative - exact_derivative) /
                         sum_of_absolute_derivative_values);
          }

        deallog << "Length ratio " << length_ratio << ", degree " << degree
                << ", coefficients:";
        for (const double coefficient : coefficients)
          deallog << " " << coefficient;
        deallog << std::endl;
        deallog << "  Polynomials of degree " << degree
                << " are extrapolated exactly: "
                << (maximal_relative_error < extrapolation_tolerance ? "true" :
                                                                       "false")
                << std::endl;
        deallog << "  Derivative coefficients:";
        for (const double coefficient : derivative_coefficients)
          deallog << " " << coefficient;
        deallog << std::endl;
        deallog << "  Derivatives of the polynomials of degree " << degree
                << " are exact: "
                << (maximal_derivative_relative_error <
                        extrapolation_tolerance ?
                      "true" :
                      "false")
                << std::endl;
      }

  // A disk of radius 0.5 centered at the origin, and a degree of freedom at a
  // distance 0.5 from its surface.
  IBParticle<2> particle;
  particle.initialize_all();
  particle.initialize_shape("sphere", std::vector<double>{0.5});
  const Point<2> dof_point(0.6, 0.8);

  for (unsigned int degree = 1; degree <= 2; ++degree)
    {
      const auto [point, interpolation_points] =
        stencil.support_points_for_interpolation(degree,
                                                 4.,
                                                 particle,
                                                 dof_point);
      deallog << "Degree " << degree << ", stencil point: " << point[0] << " "
              << point[1] << ", interpolation points:";
      for (const Point<2> &interpolation_point : interpolation_points)
        deallog << " (" << interpolation_point[0] << " "
                << interpolation_point[1] << ")";
      deallog << std::endl;
    }

  const Point<2> detection_point =
    stencil.point_for_cell_detection(particle, dof_point);
  deallog << "Point for the cell detection: " << detection_point[0] << " "
          << detection_point[1] << std::endl;

  // Derivative of f(x, y) = x^2 + y along the normal pointing out of the
  // particle, at the point of the surface closest to the degree of freedom,
  // (0.3, 0.4). The degree of freedom is at a distance d = 0.5 from the
  // surface, and the derivative is -length_ratio / d times the derivative of
  // the stencil polynomial (see IBStencil::derivative_coefficients()).
  const auto quadratic_function = [](const Point<2> &point) {
    return point[0] * point[0] + point[1];
  };
  const double       length_ratio = 4.;
  const unsigned int degree       = 2;
  const auto [point, interpolation_points] =
    stencil.support_points_for_interpolation(degree,
                                             length_ratio,
                                             particle,
                                             dof_point);
  const std::vector<double> derivative_coefficients =
    stencil.derivative_coefficients(degree, length_ratio);
  double stencil_derivative =
    derivative_coefficients[0] * quadratic_function(point);
  for (unsigned int k = 1; k <= degree; ++k)
    stencil_derivative += derivative_coefficients[k] *
                          quadratic_function(interpolation_points[k - 1]);
  const double distance_to_surface = 0.5;
  deallog << "Degree " << degree
          << ", normal derivative of x^2 + y on the surface: "
          << -length_ratio / distance_to_surface * stencil_derivative
          << " (exact: 1.16)" << std::endl;
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
