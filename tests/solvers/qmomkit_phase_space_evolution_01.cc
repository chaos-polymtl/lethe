// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the coupling between Lethe and QMOMKit on a complete QMOM
 * case, which goes through the whole QMOMKit pipeline: the quadrature of a
 * distribution is reconstructed from its moments (Wheeler algorithm), evolved
 * in phase space (node advection and moment source, combined by Strang or Lie
 * splitting), and the resulting moments are compared against a closed-form
 * solution.
 *
 * The velocity is dxi/dt = xi, so each abscissa grows as xi(0)*exp(t), the
 * weights are left untouched by the advection, and the raw moment of order k
 * follows M_k(t) = M_k(0)*exp(k*t). Adding a source that injects 0.1 into
 * every moment per unit time gives M_0(t) = M_0(0) + 0.1*t and, for k >= 1,
 * M_k(t) = M_k(0)*exp(k*t) + 0.1*(exp(k*t) - 1)/k. The initial moments are
 * those of three chosen nodes, so the initial state is exact as well.
 *
 * Three cases are verified:
 * - the advection alone moves the abscissas to xi(0)*exp(t), leaves the
 *   weights unchanged and the rebuilt moments match the exact ones;
 * - with the moment source, the error decreases as the time step is refined,
 *   at second order for the Strang splitting and first order for the two Lie
 *   splittings;
 * - a source written through the node visitor of the source context, which
 *   grows every node at a constant rate, translates the nodes exactly.
 *
 * This test is adapted from the QMOMKit verification test
 * tests/verification/phase_space_evolution/qmom_linear_advection_manufactured.cpp
 * (QMOMKit v1.0.0, Apache-2.0, Copyright 2026 The QMOMKit developers).
 */

// Deal.II
#include <deal.II/base/exceptions.h>

// QMOMKit
#include <qmomkit/moment_reconstruction.hpp>
#include <qmomkit/moment_set.hpp>
#include <qmomkit/moment_source.hpp>
#include <qmomkit/node_advection.hpp>
#include <qmomkit/phase_space_evolution.hpp>
#include <qmomkit/qmom.hpp>
#include <qmomkit/quadrature.hpp>

// Tests
#include <../tests/tests.h>

// Standard
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace
{
  /// Number of nodes of the quadrature
  constexpr std::size_t node_count = 3;

  /// Time at which the solution is compared against the exact one
  constexpr double final_time = 0.1;

  /// Rate at which the source injects into every moment
  constexpr double injection_rate = 0.1;

  /**
   * @brief Throws if a verification check fails.
   *
   * @param[in] condition Result of the check.
   * @param[in] description Description of the check, reported on failure.
   */
  void
  verify(const bool condition, const std::string &description)
  {
    AssertThrow(condition,
                ExcMessage("QMOMKit verification failed: " + description));
  }

  /**
   * @brief Returns value^exponent for an integer exponent.
   *
   * @param[in] value Base.
   * @param[in] exponent Exponent.
   */
  double
  integer_power(const double value, const std::size_t exponent)
  {
    double result = 1.0;
    for (std::size_t i = 0; i < exponent; ++i)
      result *= value;
    return result;
  }

  /**
   * @brief Returns the raw moments of the three-node distribution of the test.
   *
   * @param[in] count Number of moments.
   */
  std::vector<double>
  initial_moments(const std::size_t count)
  {
    constexpr std::array<double, node_count> nodes{{0.5, 1.5, 3.0}};
    constexpr std::array<double, node_count> weights{{0.2, 0.3, 0.5}};
    std::vector<double>                      moments(count, 0.0);
    for (std::size_t order = 0; order < count; ++order)
      for (std::size_t node = 0; node < node_count; ++node)
        moments[order] += weights[node] * integer_power(nodes[node], order);
    return moments;
  }

  /**
   * @brief Returns the options of the QMOM inversion, with a fixed number of
   * nodes.
   */
  qmomkit::QmomOptions
  qmom_options()
  {
    qmomkit::QmomOptions options;
    options.adaptive = false;
    return options;
  }

  /**
   * @brief Returns the QMOM state whose quadrature is reconstructed from
   * moments.
   *
   * @param[in] moments Raw moments.
   * @param[in] nodes Number of nodes of the quadrature.
   */
  qmomkit::QmomState
  make_state(const std::vector<double> &moments, const std::size_t nodes)
  {
    qmomkit::QmomCase qmom_case;
    qmom_case.moments      = qmomkit::MomentSet(moments);
    qmom_case.node_count   = nodes;
    qmom_case.qmom_options = qmom_options();
    return qmom_case.initial_state();
  }

  /**
   * @brief Returns the QMOM state built from the initial moments.
   *
   * @param[in] count Number of moments.
   */
  qmomkit::QmomState
  initial_state(const std::size_t count)
  {
    return make_state(initial_moments(count), node_count);
  }

  /**
   * @brief Returns the phase-space velocity dxi/dt = xi.
   */
  qmomkit::PhaseSpaceVelocity
  velocity()
  {
    return [](double,
              const std::vector<double> &coordinates,
              std::vector<double>       &velocities,
              std::size_t) {
      std::copy(coordinates.begin(), coordinates.end(), velocities.begin());
    };
  }

  /**
   * @brief Returns the maximum error between two sets of moments, relative to
   * max(|expected|, 1).
   *
   * @param[in] actual Computed moments.
   * @param[in] expected Exact moments.
   */
  double
  max_relative_error(const std::vector<double> &actual,
                     const std::vector<double> &expected)
  {
    verify(actual.size() == expected.size(),
           "the number of moments matches the exact solution");
    double result = 0.0;
    for (std::size_t i = 0; i < actual.size(); ++i)
      result = std::max(result,
                        std::abs(actual[i] - expected[i]) /
                          std::max(1.0, std::abs(expected[i])));
    return result;
  }

  /**
   * @brief Returns the options of the phase-space evolution.
   *
   * @param[in] scheme Splitting scheme combining the advection and the source.
   */
  qmomkit::QmomPhaseSpaceEvolutionOptions
  evolution_options(const qmomkit::PhaseSpaceEvolutionScheme scheme)
  {
    qmomkit::QmomPhaseSpaceEvolutionOptions result;
    result.node_options.scheme = qmomkit::NodeAdvectionScheme::rk45;
    result.node_options.adaptive_options.relative_tolerance = 1.0e-12;
    result.node_options.adaptive_options.absolute_tolerance = 1.0e-14;
    result.node_options.adaptive_options.initial_step       = 0.0025;
    result.node_options.adaptive_options.max_step           = 0.0025;
    result.scheme                                           = scheme;
    return result;
  }

  /**
   * @brief Logs the abscissas and weights of a quadrature.
   *
   * @param[in] quadrature Quadrature to log.
   */
  void
  log_quadrature(const qmomkit::QmomResult &quadrature)
  {
    for (std::size_t node = 0; node < quadrature.abscissas.size(); ++node)
      deallog << "  Node " << node
              << ": abscissa = " << quadrature.abscissas[node]
              << ", weight = " << quadrature.weights[node] << std::endl;
  }

  /**
   * @brief Advects the quadrature alone and compares the abscissas, the
   * weights and the rebuilt moments against the exact solution.
   */
  void
  test_direct_phase_space_advection()
  {
    constexpr std::size_t moment_count       = 11;
    constexpr double      time_step          = 0.00025;
    constexpr std::size_t step_count         = 400;
    constexpr double      abscissa_tolerance = 2.0e-12;
    constexpr double      weight_tolerance   = 1.0e-15;
    constexpr double      moment_tolerance   = 2.0e-11;

    deallog << "Advection of the quadrature" << std::endl;

    qmomkit::QmomState        state              = initial_state(moment_count);
    const qmomkit::QmomResult initial_quadrature = state.quadrature();
    qmomkit::QmomPhaseSpaceEvolutionOptions advection_options =
      evolution_options(qmomkit::PhaseSpaceEvolutionScheme::strang);
    advection_options.node_options.scheme =
      qmomkit::NodeAdvectionScheme::ssp_rk3;
    for (std::size_t step = 0; step < step_count; ++step)
      {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
          qmomkit::advance_qmom_advection(state,
                                          static_cast<double>(step) * time_step,
                                          time_step,
                                          velocity(),
                                          advection_options);
        verify(result.ok(), "the advection step succeeds: " + result.message);
        state = result.state;
      }

    log_quadrature(state.quadrature());

    for (std::size_t node = 0; node < node_count; ++node)
      {
        verify(std::abs(state.quadrature().abscissas[node] -
                        initial_quadrature.abscissas[node] *
                          std::exp(final_time)) < abscissa_tolerance,
               "the abscissas grow as xi(0)*exp(t)");
        verify(std::abs(state.quadrature().weights[node] -
                        initial_quadrature.weights[node]) < weight_tolerance,
               "the advection leaves the weights unchanged");
      }
    deallog << "  Abscissas grow as xi(0)*exp(t): passed" << std::endl;
    deallog << "  Weights are unchanged: passed" << std::endl;

    std::vector<double> expected = initial_moments(moment_count);
    for (std::size_t order = 0; order < expected.size(); ++order)
      expected[order] *= std::exp(static_cast<double>(order) * final_time);
    const qmomkit::Quadrature1D final_quadrature{
      state.quadrature().abscissas,
      state.quadrature().weights,
      qmomkit::QuadratureStatus::success,
      ""};
    verify(max_relative_error(qmomkit::reconstruct_raw_moments(final_quadrature,
                                                               expected.size()),
                              expected) < moment_tolerance,
           "the rebuilt moments match M_k(0)*exp(k*t)");
    deallog << "  Rebuilt moments match M_k(0)*exp(k*t): passed" << std::endl;
  }

  /**
   * @brief Returns the source that injects injection_rate into every moment.
   */
  qmomkit::QmomMomentSource
  injection_source()
  {
    return [](double, const qmomkit::QmomSourceContext &context) {
      return qmomkit::QmomMomentRate{
        std::vector<double>(context.size(), injection_rate)};
    };
  }

  /**
   * @brief Returns the exact moments at the final time with the injection
   * source.
   *
   * @param[in] count Number of moments.
   */
  std::vector<double>
  exact_injection_solution(const std::size_t count)
  {
    std::vector<double> result = initial_moments(count);
    result[0] += injection_rate * final_time;
    for (std::size_t order = 1; order < result.size(); ++order)
      {
        const double k = static_cast<double>(order);
        result[order]  = result[order] * std::exp(k * final_time) +
                        injection_rate * (std::exp(k * final_time) - 1.0) / k;
      }
    return result;
  }

  /**
   * @brief Evolves the initial state with the velocity and the injection
   * source up to the final time and returns its moments.
   *
   * @param[in] count Number of moments.
   * @param[in] time_step Time step.
   * @param[in] scheme Splitting scheme combining the advection and the source.
   */
  std::vector<double>
  evolve(const std::size_t                        count,
         const double                             time_step,
         const qmomkit::PhaseSpaceEvolutionScheme scheme)
  {
    qmomkit::QmomState state = initial_state(count);
    const std::size_t  step_count =
      static_cast<std::size_t>(std::lround(final_time / time_step));
    for (std::size_t step = 0; step < step_count; ++step)
      {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
          qmomkit::advance_qmom_phase_space(state,
                                            static_cast<double>(step) *
                                              time_step,
                                            time_step,
                                            velocity(),
                                            injection_source(),
                                            evolution_options(scheme));
        verify(result.ok(),
               "the phase-space evolution step succeeds: " + result.message);
        verify(result.diagnostics.moment_source_diagnostics.accepted,
               "the moment source step is accepted");
        state = result.state;
      }
    return state.moment_values();
  }

  /**
   * @brief Evolves the quadrature with the velocity and the injection source,
   * and checks that each splitting scheme converges to the exact solution at
   * its expected order.
   */
  void
  test_phase_space_source_evolution()
  {
    constexpr std::size_t moment_count = 6;
    // Bounds of the observed orders, which are only approximations of the
    // orders of the schemes since they are measured from two finite steps
    constexpr double min_second_order = 1.7;
    constexpr double min_first_order  = 0.7;
    constexpr double max_first_order  = 1.3;

    deallog << "Phase-space evolution with a moment source" << std::endl;

    const std::vector<double> time_steps{0.02, 0.01, 0.005, 0.0025};
    const std::array<std::pair<qmomkit::PhaseSpaceEvolutionScheme, std::string>,
                     3>
      schemes{{{qmomkit::PhaseSpaceEvolutionScheme::strang, "strang"},
               {qmomkit::PhaseSpaceEvolutionScheme::advsrc, "advsrc"},
               {qmomkit::PhaseSpaceEvolutionScheme::srcadv, "srcadv"}}};
    const std::vector<double> exact = exact_injection_solution(moment_count);

    std::array<std::vector<double>, 3> errors;
    for (std::size_t scheme = 0; scheme < schemes.size(); ++scheme)
      {
        for (const double time_step : time_steps)
          errors[scheme].push_back(max_relative_error(
            evolve(moment_count, time_step, schemes[scheme].first), exact));
        for (std::size_t run = 1; run < errors[scheme].size(); ++run)
          verify(errors[scheme][run] < errors[scheme][run - 1],
                 "the error of " + schemes[scheme].second +
                   " decreases as the time step is refined");
      }
    deallog << "  Errors decrease as the time step is refined: passed"
            << std::endl;

    // Observed orders over the last two refinements
    std::array<std::vector<double>, 3> orders;
    for (std::size_t scheme = 0; scheme < schemes.size(); ++scheme)
      for (std::size_t run = errors[scheme].size() - 2;
           run < errors[scheme].size();
           ++run)
        orders[scheme].push_back(
          std::log(errors[scheme][run - 1] / errors[scheme][run]) /
          std::log(2.0));

    // The orders are rounded to one decimal, which is all the precision two
    // finite steps give, so that the output does not depend on the platform
    for (std::size_t scheme = 0; scheme < schemes.size(); ++scheme)
      {
        std::ostringstream rounded_orders;
        rounded_orders << std::fixed << std::setprecision(1);
        for (const double order : orders[scheme])
          rounded_orders << " " << order;
        deallog << "  Observed orders of " << schemes[scheme].second << ":"
                << rounded_orders.str() << std::endl;
      }

    for (const double order : orders[0])
      verify(order > min_second_order, "strang is second order");
    for (std::size_t scheme = 1; scheme < schemes.size(); ++scheme)
      for (const double order : orders[scheme])
        verify(order > min_first_order && order < max_first_order,
               schemes[scheme].second + " is first order");
    deallog << "  Strang splitting is second order: passed" << std::endl;
    deallog << "  Lie splittings are first order: passed" << std::endl;
  }

  /**
   * @brief Grows every node at a constant rate with a source written through
   * the node visitor of the source context, and checks that the nodes are
   * translated exactly.
   */
  void
  test_node_source_translation()
  {
    constexpr std::array<double, 2> nodes{{0.5, 2.0}};
    constexpr std::array<double, 2> weights{{0.4, 0.6}};
    constexpr double                growth_rate = 0.75;
    constexpr std::size_t           step_count  = 10;
    constexpr double                time_step   = final_time / step_count;
    // The steps are exact, so the tolerances only leave room for roundoff
    constexpr double moment_tolerance = 1.0e-13;
    constexpr double node_tolerance   = 1.0e-12;

    deallog << "Moment source through the node visitor" << std::endl;

    const auto moments_of = [&](const double shift) {
      std::vector<double> moments(2 * nodes.size(), 0.0);
      for (std::size_t order = 0; order < moments.size(); ++order)
        for (std::size_t node = 0; node < nodes.size(); ++node)
          moments[order] +=
            weights[node] * integer_power(nodes[node] + shift, order);
      return moments;
    };

    const qmomkit::QmomMomentSource source =
      [](double, const qmomkit::QmomSourceContext &context) {
        std::vector<double> rate(context.size(), 0.0);
        context.for_each_node(
          [&rate](const double weight, const double abscissa) {
            for (std::size_t order = 1; order < rate.size(); ++order)
              rate[order] += weight * static_cast<double>(order) * growth_rate *
                             integer_power(abscissa, order - 1);
          });
        return qmomkit::QmomMomentRate{std::move(rate)};
      };

    qmomkit::QmomState state = make_state(moments_of(0.0), nodes.size());
    for (std::size_t step = 0; step < step_count; ++step)
      {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
          qmomkit::advance_qmom_moment_source(
            state, static_cast<double>(step) * time_step, time_step, source);
        verify(result.ok(),
               "the moment source step succeeds: " + result.message);
        state = result.state;
      }

    log_quadrature(state.quadrature());

    const double shift = growth_rate * final_time;
    verify(max_relative_error(state.moment_values(), moments_of(shift)) <
             moment_tolerance,
           "the moments are those of the translated nodes");
    for (std::size_t node = 0; node < nodes.size(); ++node)
      {
        verify(std::abs(state.quadrature().abscissas[node] -
                        (nodes[node] + shift)) < node_tolerance,
               "the abscissas are translated by the growth");
        verify(std::abs(state.quadrature().weights[node] - weights[node]) <
                 node_tolerance,
               "the weights are unchanged by the growth");
      }
    deallog << "  Moments, abscissas and weights match the translated nodes: "
               "passed"
            << std::endl;
  }
} // namespace

void
test()
{
  test_direct_phase_space_advection();
  test_phase_space_source_evolution();
  test_node_source_translation();
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
