// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_field_function.h
 * @brief Vector-valued function of the fields of the quadrature method of
 * moments (QMOM) physics: the scalars and the raw moments of a particle size
 * distribution, computed with the QMOMKit library.
 */

#ifndef lethe_qmom_field_function_h
#define lethe_qmom_field_function_h

#include <core/parameters_multiphysics.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

using namespace dealii;

/**
 * @brief Exception thrown when the QMOM physics is enabled but Lethe was not
 * linked with QMOMKit.
 */
DeclExceptionMsg(
  QMOMWithoutQMOMKitError,
  "The \"qmom\" physics is enabled, but Lethe was compiled without QMOMKit. "
  "Reconfigure Lethe with -DLETHE_WITH_QMOMKIT=ON and "
  "-DQMOMKIT_DIR=/path/to/qmomkit/installation to use this physics.");

/**
 * @brief Function of all the QMOM fields, built from the functions of the
 * scalars and from the particle size distribution.
 *
 * It has one component per field, n_scalars + n_moments components in total.
 * The scalars come first, then the moments:
 *  - components 0 to n_scalars - 1 are the scalars: component c is scalar c;
 *  - components n_scalars to n_scalars + n_moments - 1 are the moments:
 *    component c is the moment M_k, with k = c - n_scalars.
 *
 * For example, with 2 scalars and 4 moments, components 0 and 1 are the
 * scalars 0 and 1, and components 2, 3, 4 and 5 are the moments M_0, M_1, M_2
 * and M_3. Without scalars, component c is the moment M_c.
 *
 * Scalar c is given by expression c of the scalars function, so the scalar
 * expressions must be given in the same order as the scalar names.
 *
 * The moments are the raw moments of the particle size distribution, scaled by
 * the number density: M_k = n E[L^k], where n is the number density in
 * particles per m^3 and L the particle size in m. They are computed
 * analytically by QMOMKit, and are all zero if the distribution type is none.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class QMOMFieldFunction : public Function<dim>
{
public:
  /**
   * @brief Constructor. Verifies that the scalars function has one component
   * per scalar, unless there are no scalars, and that QMOMKit is available if
   * the moments are computed from a particle size distribution.
   *
   * @param[in] functions Functions of the scalars and particle size
   * distribution.
   * @param[in] n_scalars Number of scalars transported by the QMOM physics.
   * @param[in] n_moments Number of moments transported by the QMOM physics.
   * @param[in] description Description of the functions (e.g. "initial
   * conditions"), used in the error message.
   */
  QMOMFieldFunction(const Parameters::QMOMFunctions<dim> &functions,
                    const unsigned int                    n_scalars,
                    const unsigned int                    n_moments,
                    const std::string                    &description);

  /**
   * @brief Evaluate one component of the function: scalar component if
   * component < n_scalars, moment M_{component - n_scalars} otherwise.
   *
   * Evaluating a moment computes all the moments at p. Use vector_value() to
   * evaluate several moments at the same point.
   *
   * @param[in] p Point at which the function is evaluated.
   * @param[in] component Index of the field, between 0 and
   * n_scalars + n_moments - 1 (see the class description).
   *
   * @return Value of the component at p.
   */
  double
  value(const Point<dim> &p, const unsigned int component) const override;

  /**
   * @brief Evaluate all the components of the function. The moments are
   * computed once for all of them.
   *
   * @param[in] p Point at which the function is evaluated.
   * @param[out] values Values of all the components at p, of size
   * n_scalars + n_moments: values[i] is scalar i for i < n_scalars and
   * values[n_scalars + k] is the moment M_k.
   */
  void
  vector_value(const Point<dim> &p, Vector<double> &values) const override;

  /**
   * @brief Set the time of the function and of the underlying functions.
   *
   * @param[in] new_time New time.
   */
  void
  set_time(const double new_time) override;

private:
  /**
   * @brief Compute the raw moments of the particle size distribution at a
   * point.
   *
   * @param[in] p Point at which the distribution is evaluated.
   *
   * @return Raw moments M_0 to M_{n_moments - 1}.
   */
  std::vector<double>
  compute_moments(const Point<dim> &p) const;

  /// Functions of the scalars and particle size distribution
  Parameters::QMOMFunctions<dim> functions;

  /// Number of scalars
  const unsigned int n_scalars;

  /// Number of moments
  const unsigned int n_moments;
};

#endif
