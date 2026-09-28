// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_filter_kernels_h
#define lethe_filter_kernels_h

#include <cmath>

/**
 * @brief Radially symmetric Gaussian filter kernel of unit mass with compact
 * support.
 *
 * The Gaussian of standard deviation \f$\sigma\f$ is truncated at the support
 * radius \f$R = c\sigma\f$, where \f$c\f$ is the cutoff, and is renormalized
 * so that its integral over the ball of radius \f$R\f$ is exactly one:
 * \f[
 *   g(r) = \frac{1}{Z} \exp\left(-\frac{r^2}{2\sigma^2}\right)
 *   \quad \text{if } r < R, \qquad g(r) = 0 \quad \text{otherwise},
 * \f]
 * with \f$Z = (2\pi\sigma^2)^{d/2} m(c)\f$, where \f$m(c)\f$ is the mass of
 * the untruncated Gaussian contained in the support:
 * \f$m(c) = 1-\exp(-c^2/2)\f$ in 2D and
 * \f$m(c) = \mathrm{erf}(c/\sqrt{2}) - \sqrt{2/\pi}\,c\,\exp(-c^2/2)\f$ in 3D.
 *
 * Renormalizing the truncated kernel guarantees a unit kernel mass for any
 * cutoff. Without it, the mass would only be 0.971 in 3D for a cutoff of 3,
 * which would bias every filtered volume fraction by the same amount.
 *
 * @tparam dim Number of spatial dimensions (2 or 3).
 */
template <int dim>
class GaussianFilterKernel
{
  static_assert(dim == 2 || dim == 3,
                "The Gaussian filter kernel is only defined in 2D and 3D.");

public:
  /**
   * @brief Constructor.
   *
   * @param[in] standard_deviation Standard deviation of the Gaussian, i.e. the
   * physical width of the filter. Must be strictly positive.
   *
   * @param[in] cutoff Truncation radius of the kernel, expressed in number of
   * standard deviations. Must be strictly positive.
   */
  GaussianFilterKernel(const double standard_deviation, const double cutoff);

  /**
   * @brief Evaluate the kernel from the squared distance to its center. The
   * squared distance is used to avoid a square root in the hot loops of the
   * filters.
   *
   * @param[in] r_squared Squared distance between the center of the kernel and
   * the evaluation point.
   *
   * @return Value of the kernel, which is zero outside of its support.
   */
  inline double
  value_from_squared_distance(const double r_squared) const
  {
    return (r_squared < cutoff_radius_squared) ?
             normalization * std::exp(-r_squared * inverse_two_sigma_squared) :
             0.;
  }

  /**
   * @brief Return the standard deviation of the Gaussian.
   *
   * @return Standard deviation.
   */
  double
  standard_deviation() const
  {
    return sigma;
  }

  /**
   * @brief Return the radius of the support of the kernel.
   *
   * @return Support radius, i.e. the cutoff times the standard deviation.
   */
  double
  support_radius() const
  {
    return cutoff_radius;
  }

  /**
   * @brief Return the squared radius of the support of the kernel.
   *
   * @return Squared support radius.
   */
  double
  support_radius_squared() const
  {
    return cutoff_radius_squared;
  }

  /**
   * @brief Return the mass of the untruncated Gaussian contained in the
   * support. This is the fraction of the Gaussian that the truncation keeps
   * before the renormalization.
   *
   * @return Retained mass fraction, between 0 and 1.
   */
  double
  retained_mass_fraction() const
  {
    return mass_fraction;
  }

private:
  /// Standard deviation of the Gaussian.
  double sigma;

  /// Radius of the support of the kernel.
  double cutoff_radius;

  /// Squared radius of the support of the kernel.
  double cutoff_radius_squared;

  /// Inverse of twice the variance, 1/(2 sigma^2).
  double inverse_two_sigma_squared;

  /// Mass of the untruncated Gaussian contained in the support.
  double mass_fraction;

  /// Normalization constant 1/Z of the truncated kernel.
  double normalization;
};

/**
 * @brief Top-hat filter kernel, i.e. the normalized indicator function of a
 * ball:
 * \f[
 *   g(r) = \frac{1}{|B_R|} \quad \text{if } r < R, \qquad g(r) = 0 \quad
 *   \text{otherwise},
 * \f]
 * where \f$|B_R|\f$ is the area (2D) or the volume (3D) of the ball of radius
 * \f$R\f$.
 *
 * @tparam dim Number of spatial dimensions (2 or 3).
 */
template <int dim>
class TopHatFilterKernel
{
  static_assert(dim == 2 || dim == 3,
                "The top-hat filter kernel is only defined in 2D and 3D.");

public:
  /**
   * @brief Constructor.
   *
   * @param[in] radius Radius of the ball supporting the kernel. Must be
   * strictly positive.
   */
  explicit TopHatFilterKernel(const double radius);

  /**
   * @brief Evaluate the kernel from the squared distance to its center.
   *
   * @param[in] r_squared Squared distance between the center of the kernel and
   * the evaluation point.
   *
   * @return Value of the kernel, which is zero outside of its support.
   */
  inline double
  value_from_squared_distance(const double r_squared) const
  {
    return (r_squared < ball_radius_squared) ? normalization : 0.;
  }

  /**
   * @brief Return the radius of the support of the kernel.
   *
   * @return Radius of the ball.
   */
  double
  support_radius() const
  {
    return ball_radius;
  }

  /**
   * @brief Return the squared radius of the support of the kernel.
   *
   * @return Squared radius of the ball.
   */
  double
  support_radius_squared() const
  {
    return ball_radius_squared;
  }

private:
  /// Radius of the ball.
  double ball_radius;

  /// Squared radius of the ball.
  double ball_radius_squared;

  /// Normalization constant, the inverse of the measure of the ball.
  double normalization;
};

#endif
