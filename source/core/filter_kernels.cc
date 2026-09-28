// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/filter_kernels.h>

#include <deal.II/base/exceptions.h>

#include <cmath>
#include <numbers>

using namespace dealii;

namespace
{
  /**
   * @brief Return the mass of the untruncated Gaussian contained in a ball of
   * radius cutoff * sigma centered on the Gaussian.
   *
   * @tparam dim Number of spatial dimensions (2 or 3).
   *
   * @param[in] cutoff Radius of the ball, expressed in number of standard
   * deviations.
   *
   * @return Mass fraction, between 0 and 1.
   */
  template <int dim>
  double
  untruncated_gaussian_mass_fraction(const double cutoff)
  {
    const double half_cutoff_squared = 0.5 * cutoff * cutoff;
    if constexpr (dim == 2)
      return -std::expm1(-half_cutoff_squared);
    else
      return std::erf(cutoff / std::numbers::sqrt2) -
             std::sqrt(2. / std::numbers::pi) * cutoff *
               std::exp(-half_cutoff_squared);
  }

  /**
   * @brief Return the measure of a ball, i.e. its area in 2D or its volume in
   * 3D.
   *
   * @tparam dim Number of spatial dimensions (2 or 3).
   *
   * @param[in] radius Radius of the ball.
   *
   * @return Measure of the ball.
   */
  template <int dim>
  double
  ball_measure(const double radius)
  {
    if constexpr (dim == 2)
      return std::numbers::pi * radius * radius;
    else
      return 4. / 3. * std::numbers::pi * radius * radius * radius;
  }
} // namespace

template <int dim>
GaussianFilterKernel<dim>::GaussianFilterKernel(const double standard_deviation,
                                                const double cutoff)
  : sigma(standard_deviation)
  , cutoff_radius(cutoff * standard_deviation)
  , cutoff_radius_squared(cutoff_radius * cutoff_radius)
  , inverse_two_sigma_squared(0.5 / (standard_deviation * standard_deviation))
  , mass_fraction(untruncated_gaussian_mass_fraction<dim>(cutoff))
  , normalization(1. / (std::pow(2. * std::numbers::pi * standard_deviation *
                                   standard_deviation,
                                 0.5 * dim) *
                        mass_fraction))
{
  AssertThrow(standard_deviation > 0.,
              ExcMessage("The standard deviation of the Gaussian filter "
                         "kernel must be strictly positive."));
  AssertThrow(cutoff > 0.,
              ExcMessage("The cutoff of the Gaussian filter kernel must be "
                         "strictly positive."));
}

template <int dim>
TopHatFilterKernel<dim>::TopHatFilterKernel(const double radius)
  : ball_radius(radius)
  , ball_radius_squared(radius * radius)
  , normalization(1. / ball_measure<dim>(radius))
{
  AssertThrow(radius > 0.,
              ExcMessage("The radius of the top-hat filter kernel must be "
                         "strictly positive."));
}

template class GaussianFilterKernel<2>;
template class GaussianFilterKernel<3>;
template class TopHatFilterKernel<2>;
template class TopHatFilterKernel<3>;
