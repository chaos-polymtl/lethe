// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom_source_terms.h
 * @brief Mechanisms that evolve the particle size distribution of the QMOM
 * physics in phase space: growth, nucleation and coalescence.
 *
 * Only the growth is implemented so far. Each growth model is a fixed
 * function that returns the phase-space velocity of the quadrature nodes,
 * with the signature of qmomkit::PhaseSpaceVelocity, so that it is passed to
 * QMOMKit as is. The function model, which parses an expression, is a generic
 * option: the fixed models are preferred, since they are evaluated without
 * parsing.
 */

#ifndef lethe_qmom_source_terms_h
#define lethe_qmom_source_terms_h

#include <core/parameters_multiphysics.h>

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

/**
 * @brief Phase-space velocity of the quadrature nodes, with the signature of
 * qmomkit::PhaseSpaceVelocity: time, coordinates of the nodes, velocities of
 * the nodes (to fill, already sized) and number of phase-space dimensions.
 * The phase space of the QMOM physics has a single dimension, the particle
 * size xi in m, so the coordinates are the abscissas of the nodes and the
 * velocities are their growth rates, in m/s.
 */
using QMOMPhaseSpaceVelocity = std::function<void(const double,
                                                  const std::vector<double> &,
                                                  std::vector<double> &,
                                                  const std::size_t)>;

/**
 * @brief Return the phase-space velocity without growth: G = 0.
 *
 * @return Phase-space velocity that sets the growth rate of every node to 0.
 */
QMOMPhaseSpaceVelocity
no_growth_velocity();

/**
 * @brief Return the phase-space velocity of a constant growth: G = a.
 *
 * @param[in] growth_coefficient Growth rate a, in m/s.
 *
 * @return Phase-space velocity that sets the growth rate of every node to a.
 */
QMOMPhaseSpaceVelocity
constant_growth_velocity(const double growth_coefficient);

/**
 * @brief Return the phase-space velocity of a linear growth: G = a xi.
 *
 * @param[in] growth_coefficient Coefficient a, in 1/s.
 *
 * @return Phase-space velocity that sets the growth rate of every node to a
 * times its abscissa.
 */
QMOMPhaseSpaceVelocity
linear_growth_velocity(const double growth_coefficient);

/**
 * @brief Return the phase-space velocity given by an expression of the
 * particle size xi and of the time t. The expression is parsed once, but
 * evaluated by the parser for every node, which is slower than the fixed
 * models.
 *
 * @param[in] function_expression Expression of the growth rate, in m/s, in
 * terms of xi, in m, and t, in s.
 *
 * @return Phase-space velocity that sets the growth rate of every node to
 * the expression evaluated at its abscissa.
 */
QMOMPhaseSpaceVelocity
function_growth_velocity(const std::string &function_expression);

/**
 * @brief Phase-space mechanisms of the QMOM physics, built from the
 * parameters of the "QMOM" subsection.
 *
 * The growth velocity is selected once, at construction, among the fixed
 * functions above, according to the model of the "growth" subsection.
 */
class QMOMSourceTerms
{
public:
  /**
   * @brief Constructor. Selects the growth velocity of the growth model, or
   * the velocity without growth if the growth is disabled.
   *
   * @param[in] qmom_parameters Parameters of the "QMOM" subsection.
   */
  explicit QMOMSourceTerms(const Parameters::QMOM &qmom_parameters);

  /**
   * @brief Indicate whether the particles grow in phase space.
   *
   * @return True if the growth is enabled.
   */
  bool
  growth_is_enabled() const
  {
    return enable_growth;
  }

  /**
   * @brief Get the phase-space velocity of the quadrature nodes, which is
   * zero if the growth is disabled.
   *
   * @return Reference to the phase-space velocity.
   */
  const QMOMPhaseSpaceVelocity &
  get_growth_velocity() const
  {
    return growth_velocity;
  }

private:
  /// Enable the growth of the particles
  const bool enable_growth;

  /// Phase-space velocity of the quadrature nodes
  QMOMPhaseSpaceVelocity growth_velocity;
};

#endif
