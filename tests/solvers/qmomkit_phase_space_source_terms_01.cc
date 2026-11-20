// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the models of the phase-space mechanisms of the QMOM physics
 * (growth, nucleation and coalescence), through the QMOMKit callbacks built by
 * QMOMSourceTerms, against closed-form solutions.
 *
 * A distribution of three nodes is evolved over a time T = 0.1 with each
 * mechanism, and its moments are compared with the exact ones:
 * - constant growth G = a: the nodes are translated by a T;
 * - linear growth G = a xi, and its function equivalent: M_k = M_k(0) e^(k a
 * T);
 * - constant nucleation J = a of nuclei whose sizes follow a log-normal
 *   distribution of moments E[xi^k], and a function nucleation J = a (1 + t):
 *   M_k = M_k(0) + E[xi^k] integral of J;
 * - linear growth with constant nucleation: M_0 = M_0(0) + J T and, for k >= 1,
 *   M_k = M_k(0) e^(k a T) + J E[xi^k] (e^(k a T) - 1) / (k a);
 * - constant coalescence kernel beta, and its function equivalent:
 *   M_0 = M_0(0) / (1 + beta M_0(0) T / 2), while M_3, proportional to the
 *   volume of the particles, is conserved.
 */

// Lethe
#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

// Deal.II
#include <deal.II/base/exceptions.h>

// QMOMKit
#include <qmomkit/moment_set.hpp>
#include <qmomkit/phase_space_evolution.hpp>
#include <qmomkit/qmom.hpp>

// Tests
#include <../tests/tests.h>

// Standard
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace
{
  /// Number of nodes of the quadrature
  constexpr std::size_t n_nodes = 3;

  /// Number of moments, which define the three nodes
  constexpr std::size_t n_moments = 2 * n_nodes;

  /// Final time of the evolution
  constexpr double final_time = 0.1;

  /// Number of time steps to reach the final time
  constexpr unsigned int n_time_steps = 100;

  /// Abscissas and weights of the initial distribution
  constexpr std::array<double, n_nodes> initial_abscissas{{0.5, 1.5, 3.0}};
  constexpr std::array<double, n_nodes> initial_weights{{0.2, 0.3, 0.5}};

  /**
   * @brief Throw if a verification fails.
   *
   * @param[in] condition Result of the verification.
   * @param[in] description Description of the verification.
   */
  void
  verify(const bool condition, const std::string &description)
  {
    AssertThrow(condition, ExcMessage("Verification failed: " + description));
  }

  /**
   * @brief Return the moments of the initial nodes, each translated by a
   * distance.
   *
   * @param[in] translation Distance by which every node is translated.
   *
   * @return Moments M_0 to M_{n_moments - 1}.
   */
  std::vector<double>
  translated_moments(const double translation)
  {
    std::vector<double> moments(n_moments, 0.);
    for (std::size_t node = 0; node < n_nodes; ++node)
      for (std::size_t order = 0; order < n_moments; ++order)
        moments[order] +=
          initial_weights[node] *
          std::pow(initial_abscissas[node] + translation, order);
    return moments;
  }

  /**
   * @brief Return the maximum error between moments, relative to each exact
   * moment.
   *
   * @param[in] moments Computed moments.
   * @param[in] exact_moments Exact moments.
   * @param[in] orders Orders of the moments that are compared.
   *
   * @return Maximum relative error.
   */
  double
  max_relative_error(const std::vector<double>       &moments,
                     const std::vector<double>       &exact_moments,
                     const std::vector<unsigned int> &orders)
  {
    double error = 0.;
    for (const unsigned int order : orders)
      error = std::max(error,
                       std::abs(moments[order] - exact_moments[order]) /
                         std::abs(exact_moments[order]));
    return error;
  }

  /**
   * @brief Evolve the initial distribution up to the final time with the
   * mechanisms of the QMOM parameters.
   *
   * @param[in] qmom_parameters Parameters of the QMOM subsection.
   *
   * @return Moments at the final time.
   */
  std::vector<double>
  evolve(const Parameters::QMOM &qmom_parameters)
  {
    QMOMSourceTerms source_terms(qmom_parameters);
    // There are no scalars
    source_terms.set_local_scalars({});

    qmomkit::QmomCase qmom_case;
    qmom_case.moments    = qmomkit::MomentSet(translated_moments(0.));
    qmom_case.node_count = n_nodes;
    qmom_case.qmom_options.adaptive = false;
    qmomkit::QmomState state        = qmom_case.initial_state();

    qmomkit::QmomPhaseSpaceEvolutionOptions evolution_options;
    evolution_options.node_options.adaptive_options.relative_tolerance = 1e-12;
    evolution_options.node_options.adaptive_options.absolute_tolerance = 1e-14;

    const double time_step = final_time / n_time_steps;
    for (unsigned int step = 0; step < n_time_steps; ++step)
      {
        const qmomkit::QmomPhaseSpaceEvolutionResult result =
          qmomkit::advance_qmom_phase_space(
            state,
            step * time_step,
            time_step,
            source_terms.get_phase_space_velocity(),
            source_terms.get_moment_source(),
            source_terms.get_boundary_contribution(),
            evolution_options);
        verify(result.ok(), "the phase-space step succeeds: " + result.message);
        state = result.state;
      }
    return state.moment_values();
  }

  /**
   * @brief Evolve the initial distribution and compare its moments with the
   * exact ones.
   *
   * @param[in] description Description of the case.
   * @param[in] qmom_parameters Parameters of the QMOM subsection.
   * @param[in] exact_moments Exact moments at the final time.
   * @param[in] orders Orders of the moments that are compared.
   * @param[in] tolerance Maximum relative error.
   */
  void
  test_case(const std::string               &description,
            const Parameters::QMOM          &qmom_parameters,
            const std::vector<double>       &exact_moments,
            const std::vector<unsigned int> &orders,
            const double                     tolerance)
  {
    const std::vector<double> moments = evolve(qmom_parameters);
    verify(max_relative_error(moments, exact_moments, orders) < tolerance,
           description + " matches the exact moments");
    deallog << description << ": passed" << std::endl;
  }

  /// All the moment orders
  const std::vector<unsigned int> all_orders{0, 1, 2, 3, 4, 5};

  void
  test_growth()
  {
    constexpr double growth_coefficient = 0.5;

    Parameters::QMOM constant_growth;
    constant_growth.enable_growth = true;
    constant_growth.growth.model  = Parameters::QMOMGrowthModel::constant;
    constant_growth.growth.growth_coefficient = growth_coefficient;
    test_case("Constant growth",
              constant_growth,
              translated_moments(growth_coefficient * final_time),
              all_orders,
              1e-10);

    std::vector<double> exponential_moments = translated_moments(0.);
    for (unsigned int order = 0; order < n_moments; ++order)
      exponential_moments[order] *=
        std::exp(order * growth_coefficient * final_time);

    Parameters::QMOM linear_growth;
    linear_growth.enable_growth = true;
    linear_growth.growth.model  = Parameters::QMOMGrowthModel::linear;
    linear_growth.growth.growth_coefficient = growth_coefficient;
    test_case(
      "Linear growth", linear_growth, exponential_moments, all_orders, 1e-10);

    Parameters::QMOM function_growth;
    function_growth.enable_growth = true;
    function_growth.growth.model  = Parameters::QMOMGrowthModel::function;
    function_growth.growth.function_expression = "0.5 * xi";
    test_case("Function growth",
              function_growth,
              exponential_moments,
              all_orders,
              1e-10);
  }

  void
  test_nucleation()
  {
    constexpr double nucleation_coefficient    = 2.;
    constexpr double nuclei_mean_size          = 0.8;
    constexpr double nuclei_standard_deviation = 0.2;

    // Moments E[xi^k] of the log-normal distribution of the nuclei
    const double nuclei_log_variance =
      std::log(1. + std::pow(nuclei_standard_deviation / nuclei_mean_size, 2));
    const double nuclei_log_mean =
      std::log(nuclei_mean_size) - 0.5 * nuclei_log_variance;
    std::vector<double> nuclei_moments(n_moments);
    for (unsigned int order = 0; order < n_moments; ++order)
      nuclei_moments[order] = std::exp(
        order * nuclei_log_mean + 0.5 * order * order * nuclei_log_variance);

    Parameters::QMOM constant_nucleation;
    constant_nucleation.enable_nucleation = true;
    constant_nucleation.nucleation.model =
      Parameters::QMOMNucleationModel::constant;
    constant_nucleation.nucleation.nucleation_coefficient =
      nucleation_coefficient;
    constant_nucleation.nucleation.nuclei_mean_size = nuclei_mean_size;
    constant_nucleation.nucleation.nuclei_standard_deviation =
      nuclei_standard_deviation;

    std::vector<double> constant_moments = translated_moments(0.);
    for (unsigned int order = 0; order < n_moments; ++order)
      constant_moments[order] +=
        nucleation_coefficient * nuclei_moments[order] * final_time;
    test_case("Constant nucleation",
              constant_nucleation,
              constant_moments,
              all_orders,
              1e-10);

    // J = a (1 + t), whose integral over [0, T] is a (T + T^2 / 2)
    Parameters::QMOM function_nucleation;
    function_nucleation.enable_nucleation = true;
    function_nucleation.nucleation.model =
      Parameters::QMOMNucleationModel::function;
    function_nucleation.nucleation.function_expression = "2 * (1 + t)";
    function_nucleation.nucleation.nuclei_mean_size    = nuclei_mean_size;
    function_nucleation.nucleation.nuclei_standard_deviation =
      nuclei_standard_deviation;

    std::vector<double> function_moments = translated_moments(0.);
    for (unsigned int order = 0; order < n_moments; ++order)
      function_moments[order] += nucleation_coefficient *
                                 nuclei_moments[order] *
                                 (final_time + 0.5 * final_time * final_time);
    test_case("Function nucleation",
              function_nucleation,
              function_moments,
              all_orders,
              1e-10);

    // Linear growth with a constant nucleation. The splitting between the two
    // mechanisms inside the phase-space step is second order in time.
    constexpr double growth_coefficient    = 0.5;
    Parameters::QMOM growth_and_nucleation = constant_nucleation;
    growth_and_nucleation.enable_growth    = true;
    growth_and_nucleation.growth.model = Parameters::QMOMGrowthModel::linear;
    growth_and_nucleation.growth.growth_coefficient = growth_coefficient;

    std::vector<double> combined_moments = translated_moments(0.);
    combined_moments[0] += nucleation_coefficient * final_time;
    for (unsigned int order = 1; order < n_moments; ++order)
      {
        const double growth_factor =
          std::exp(order * growth_coefficient * final_time);
        combined_moments[order] =
          combined_moments[order] * growth_factor +
          nucleation_coefficient * nuclei_moments[order] *
            (growth_factor - 1.) / (order * growth_coefficient);
      }
    test_case("Linear growth with constant nucleation",
              growth_and_nucleation,
              combined_moments,
              all_orders,
              1e-6);
  }

  void
  test_coalescence()
  {
    constexpr double coalescence_coefficient = 0.5;

    // Only M_0 and M_3 have a closed form with a constant kernel
    std::vector<double> exact_moments = translated_moments(0.);
    exact_moments[0] /=
      1. + 0.5 * coalescence_coefficient * exact_moments[0] * final_time;
    const std::vector<unsigned int> closed_form_orders{0, 3};

    Parameters::QMOM constant_coalescence;
    constant_coalescence.enable_coalescence = true;
    constant_coalescence.coalescence.model =
      Parameters::QMOMCoalescenceModel::constant;
    constant_coalescence.coalescence.coalescence_coefficient =
      coalescence_coefficient;
    test_case("Constant coalescence",
              constant_coalescence,
              exact_moments,
              closed_form_orders,
              1e-8);

    Parameters::QMOM function_coalescence;
    function_coalescence.enable_coalescence = true;
    function_coalescence.coalescence.model =
      Parameters::QMOMCoalescenceModel::function;
    function_coalescence.coalescence.function_expression =
      "0.5 + 0 * (xi_1 + xi_2)";
    test_case("Function coalescence",
              function_coalescence,
              exact_moments,
              closed_form_orders,
              1e-8);
  }
} // namespace

int
main()
{
  try
    {
      initlog();
      test_growth();
      test_nucleation();
      test_coalescence();
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
