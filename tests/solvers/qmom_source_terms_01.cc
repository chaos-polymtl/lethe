// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @brief Tests the growth, nucleation and coalescence models of the QMOM
 * physics against their exact formulas:
 * - growth: G = a, G = a xi and G = 0.5 xi + t;
 * - nucleation: dM_k/dt = J E[xi^k], with log-normal nuclei;
 * - coalescence of a single node: dM_k/dt = 1/2 w^2 beta xi^k (2^(k/3) - 2),
 *   and, with three nodes, dM_0/dt = -beta M_0^2 / 2 and dM_3/dt = 0.
 */

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

#include <deal.II/base/utilities.h>

#include <../tests/tests.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

namespace
{
  /// Number of checked moments
  constexpr unsigned int n_moments = 4;

  /// Mean size and standard deviation of the nuclei
  constexpr double nuclei_mean_size          = 2.;
  constexpr double nuclei_standard_deviation = 0.5;

  /// Round-off tolerance on the rates
  constexpr double tolerance = 1e-12;

  /**
   * @brief Print the growth rates of three nodes and check them.
   *
   * @param[in] label Printed label.
   * @param[in] growth_parameters Growth parameters.
   * @param[in] time Time.
   * @param[in] exact_growth_rate Exact formula G(xi, t).
   */
  void
  check_growth_model(
    const std::string                           &label,
    const Parameters::QMOMGrowth                &growth_parameters,
    const double                                 time,
    const std::function<double(double, double)> &exact_growth_rate)
  {
    const std::vector<double> abscissas{0.5, 1.5, 3.0};
    const std::vector<double> local_scalars;
    std::vector<double>       growth_rates(abscissas.size());

    const auto growth_model =
      QMOMGrowthModelBase::model_cast(growth_parameters);
    growth_model->compute_growth_rates(time,
                                       local_scalars,
                                       abscissas,
                                       growth_rates);

    deallog << label << ", t = " << time << std::endl;
    bool matches_exact_rate = true;
    for (unsigned int node = 0; node < abscissas.size(); ++node)
      {
        deallog << "  xi = " << abscissas[node]
                << ": G = " << growth_rates[node] << std::endl;
        matches_exact_rate &=
          std::abs(growth_rates[node] -
                   exact_growth_rate(abscissas[node], time)) < tolerance;
      }
    deallog << "  matches the exact formula: "
            << (matches_exact_rate ? "true" : "false") << std::endl;
  }

  /**
   * @brief Compute the moment E[xi^k] of the log-normal nuclei.
   *
   * @param[in] order Order k.
   *
   * @return E[xi^k].
   */
  double
  nuclei_moment(const unsigned int order)
  {
    const double variance_ratio =
      1. +
      Utilities::fixed_power<2>(nuclei_standard_deviation / nuclei_mean_size);
    return std::pow(nuclei_mean_size, order) *
           std::pow(variance_ratio, 0.5 * order * (order - 1.));
  }

  /**
   * @brief Print the nucleation rates of the moments and check them.
   *
   * @param[in] label Printed label.
   * @param[in] moment_rates Rates to check.
   * @param[in] nucleation_rate Nucleation rate J.
   */
  void
  check_nucleation_moment_rates(const std::string         &label,
                                const std::vector<double> &moment_rates,
                                const double               nucleation_rate)
  {
    deallog << label << std::endl;
    bool matches_exact_rates = true;
    for (unsigned int order = 0; order < n_moments; ++order)
      {
        deallog << "  dM_" << order << "/dt = " << moment_rates[order]
                << std::endl;
        matches_exact_rates &=
          std::abs(moment_rates[order] -
                   nucleation_rate * nuclei_moment(order)) <
          tolerance * std::max(1., std::abs(moment_rates[order]));
      }
    deallog << "  matches the exact formula: "
            << (matches_exact_rates ? "true" : "false") << std::endl;
  }

  /**
   * @brief Print the coalescence rates of a single node and check them.
   *
   * @param[in] label Printed label.
   * @param[in] moment_rates Rates to check.
   * @param[in] weight Weight w.
   * @param[in] size Size xi.
   * @param[in] kernel Kernel beta(xi, xi).
   */
  void
  check_single_node_coalescence_rates(const std::string         &label,
                                      const std::vector<double> &moment_rates,
                                      const double               weight,
                                      const double               size,
                                      const double               kernel)
  {
    deallog << label << std::endl;
    bool matches_exact_rates = true;
    for (unsigned int order = 0; order < n_moments; ++order)
      {
        const double exact_rate = 0.5 * weight * weight * kernel *
                                  std::pow(size, order) *
                                  (std::pow(2., order / 3.) - 2.);
        // Round-off is printed as zero
        deallog << "  dM_" << order << "/dt = "
                << (std::abs(moment_rates[order]) < tolerance ?
                      0. :
                      moment_rates[order])
                << std::endl;
        matches_exact_rates &= std::abs(moment_rates[order] - exact_rate) <
                               tolerance * std::max(1., std::abs(exact_rate));
      }
    deallog << "  matches the exact formula: "
            << (matches_exact_rates ? "true" : "false") << std::endl;
  }

  /**
   * @brief Test the growth models.
   */
  void
  test_growth_models()
  {
    Parameters::QMOMGrowth growth_parameters;

    growth_parameters.model = Parameters::QMOMGrowthModel::constant;
    growth_parameters.growth_coefficient = 2.;
    check_growth_model("Constant growth, a = 2",
                       growth_parameters,
                       0.,
                       [](const double, const double) { return 2.; });

    growth_parameters.model              = Parameters::QMOMGrowthModel::linear;
    growth_parameters.growth_coefficient = 0.5;
    check_growth_model("Linear growth, a = 0.5",
                       growth_parameters,
                       0.,
                       [](const double size, const double) {
                         return 0.5 * size;
                       });

    growth_parameters.model = Parameters::QMOMGrowthModel::function;
    growth_parameters.function_expression = "0.5 * xi + t";
    check_growth_model("Function growth, G = 0.5 xi + t",
                       growth_parameters,
                       2.,
                       [](const double size, const double time) {
                         return 0.5 * size + time;
                       });
  }

  /**
   * @brief Test the nucleation models and QMOMSourceTerms.
   */
  void
  test_nucleation_models()
  {
    const std::vector<double> local_scalars;
    std::vector<double>       moment_rates(n_moments);

    Parameters::QMOMNucleation nucleation_parameters;
    nucleation_parameters.nuclei_mean_size          = nuclei_mean_size;
    nucleation_parameters.nuclei_standard_deviation = nuclei_standard_deviation;

    nucleation_parameters.model = Parameters::QMOMNucleationModel::constant;
    nucleation_parameters.nucleation_coefficient = 3.;
    QMOMNucleationModelBase::model_cast(nucleation_parameters)
      ->compute_moment_rates(0., local_scalars, moment_rates);
    check_nucleation_moment_rates("Constant nucleation, J = 3",
                                  moment_rates,
                                  3.);

    nucleation_parameters.model = Parameters::QMOMNucleationModel::function;
    nucleation_parameters.function_expression = "3 * (1 + t)";
    QMOMNucleationModelBase::model_cast(nucleation_parameters)
      ->compute_moment_rates(0.5, local_scalars, moment_rates);
    check_nucleation_moment_rates("Function nucleation, J = 3 (1 + t), t = 0.5",
                                  moment_rates,
                                  4.5);

    // QMOMSourceTerms only builds the model of an enabled mechanism
    Parameters::QMOM qmom_parameters;
    qmom_parameters.nucleation = nucleation_parameters;
    qmom_parameters.nucleation.model =
      Parameters::QMOMNucleationModel::constant;

    qmom_parameters.enable_nucleation = false;
    {
      const QMOMSourceTerms source_terms(qmom_parameters);
      source_terms.compute_nucleation_moment_rates(0., moment_rates);
      check_nucleation_moment_rates(
        "QMOMSourceTerms, nucleation disabled, nucleation_is_enabled() = " +
          std::string(source_terms.nucleation_is_enabled() ? "true" : "false"),
        moment_rates,
        0.);
    }

    qmom_parameters.enable_nucleation = true;
    {
      const QMOMSourceTerms source_terms(qmom_parameters);
      source_terms.compute_nucleation_moment_rates(0., moment_rates);
      check_nucleation_moment_rates(
        "QMOMSourceTerms, constant nucleation J = 3, nucleation_is_enabled() = " +
          std::string(source_terms.nucleation_is_enabled() ? "true" : "false"),
        moment_rates,
        3.);
    }
  }

  /**
   * @brief Test the coalescence models.
   */
  void
  test_coalescence_models()
  {
    const std::vector<double> local_scalars;
    std::vector<double>       moment_rates(n_moments);

    Parameters::QMOMCoalescence coalescence_parameters;
    coalescence_parameters.model = Parameters::QMOMCoalescenceModel::constant;
    coalescence_parameters.coalescence_coefficient = 3.;
    const auto constant_model =
      QMOMCoalescenceModelBase::model_cast(coalescence_parameters);

    coalescence_parameters.model = Parameters::QMOMCoalescenceModel::function;
    coalescence_parameters.function_expression = "(xi_1 + xi_2) * t";
    const auto function_model =
      QMOMCoalescenceModelBase::model_cast(coalescence_parameters);

    // Single node of weight 2 and size 1
    constant_model->compute_moment_rates(
      0., local_scalars, {2.}, {1.}, moment_rates);
    check_single_node_coalescence_rates(
      "Constant kernel, beta = 3, single node", moment_rates, 2., 1., 3.);

    function_model->compute_moment_rates(
      1.5, local_scalars, {2.}, {1.}, moment_rates);
    check_single_node_coalescence_rates(
      "Function kernel, beta = (xi_1 + xi_2) t, t = 1.5, single node",
      moment_rates,
      2.,
      1.,
      3.);

    constant_model->compute_moment_rates(
      0., local_scalars, {2., 5.}, {1., 0.}, moment_rates);
    check_single_node_coalescence_rates(
      "Constant kernel, beta = 3, with a second node of size zero",
      moment_rates,
      2.,
      1.,
      3.);

    // Three nodes
    const std::vector<double> weights{0.2, 0.3, 0.5};
    const std::vector<double> abscissas{0.5, 1.5, 3.0};
    constant_model->compute_moment_rates(
      0., local_scalars, weights, abscissas, moment_rates);

    const double zeroth_moment = std::reduce(weights.begin(), weights.end());
    const double third_moment =
      std::transform_reduce(weights.begin(),
                            weights.end(),
                            abscissas.begin(),
                            0.,
                            std::plus<>(),
                            [](const double weight, const double size) {
                              return weight * std::pow(size, 3);
                            });
    deallog << "Constant kernel, beta = 3, three nodes" << std::endl;
    deallog << "  dM_0/dt = " << moment_rates[0] << ", -beta M_0^2 / 2 = "
            << -0.5 * 3. * zeroth_moment * zeroth_moment << std::endl;
    deallog << "  M_3 is conserved: "
            << (std::abs(moment_rates[3]) < tolerance * third_moment ? "true" :
                                                                       "false")
            << std::endl;

    coalescence_parameters.function_expression = "3";
    std::vector<double> function_moment_rates(n_moments);
    QMOMCoalescenceModelBase::model_cast(coalescence_parameters)
      ->compute_moment_rates(
        0., local_scalars, weights, abscissas, function_moment_rates);
    bool function_matches_constant = true;
    for (unsigned int order = 0; order < n_moments; ++order)
      function_matches_constant &=
        std::abs(function_moment_rates[order] - moment_rates[order]) <
        tolerance * std::max(1., std::abs(moment_rates[order]));
    deallog << "  function kernel beta = 3 gives the same rates: "
            << (function_matches_constant ? "true" : "false") << std::endl;
  }
} // namespace

int
main()
{
  try
    {
      initlog();
      test_growth_models();
      test_nucleation_models();
      test_coalescence_models();
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
