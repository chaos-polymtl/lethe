// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace
{
  /// Constants available in the expressions of the function models
  const std::map<std::string, double> function_constants{{"pi", numbers::PI}};

} // namespace

std::unique_ptr<QMOMGrowthModelBase>
QMOMGrowthModelBase::model_cast(const Parameters::QMOMGrowth &growth_parameters)
{
  switch (growth_parameters.model)
    {
      case Parameters::QMOMGrowthModel::constant:
        return std::make_unique<QMOMConstantGrowth>(
          growth_parameters.growth_coefficient);
      case Parameters::QMOMGrowthModel::linear:
        return std::make_unique<QMOMLinearGrowth>(
          growth_parameters.growth_coefficient);
      case Parameters::QMOMGrowthModel::function:
        return std::make_unique<QMOMFunctionGrowth>(
          growth_parameters.function_expression);
    }
  AssertThrow(false, ExcInternalError());
  return nullptr;
}

QMOMConstantGrowth::QMOMConstantGrowth(const double growth_coefficient)
  : growth_coefficient(growth_coefficient)
{}

void
QMOMConstantGrowth::compute_growth_rates(
  const double /*time*/,
  const std::vector<double> & /*local_scalars*/,
  const std::vector<double> &abscissas,
  std::vector<double>       &growth_rates) const
{
  AssertDimension(growth_rates.size(), abscissas.size());
  std::ranges::fill(growth_rates, growth_coefficient);
}

QMOMLinearGrowth::QMOMLinearGrowth(const double growth_coefficient)
  : growth_coefficient(growth_coefficient)
{}

void
QMOMLinearGrowth::compute_growth_rates(
  const double /*time*/,
  const std::vector<double> & /*local_scalars*/,
  const std::vector<double> &abscissas,
  std::vector<double>       &growth_rates) const
{
  AssertDimension(growth_rates.size(), abscissas.size());
  for (unsigned int node = 0; node < abscissas.size(); ++node)
    growth_rates[node] = growth_coefficient * abscissas[node];
}

QMOMFunctionGrowth::QMOMFunctionGrowth(const std::string &function_expression)
  : growth_function(1)
{
  // The particle size xi takes the place of the spatial coordinate
  growth_function.initialize("xi,t",
                             {function_expression},
                             function_constants,
                             true);
}

void
QMOMFunctionGrowth::compute_growth_rates(
  const double time,
  const std::vector<double> & /*local_scalars*/,
  const std::vector<double> &abscissas,
  std::vector<double>       &growth_rates) const
{
  AssertDimension(growth_rates.size(), abscissas.size());
  growth_function.set_time(time);
  for (unsigned int node = 0; node < abscissas.size(); ++node)
    growth_rates[node] = growth_function.value(Point<1>(abscissas[node]));
}

std::unique_ptr<QMOMNucleationModelBase>
QMOMNucleationModelBase::model_cast(
  const Parameters::QMOMNucleation &nucleation_parameters)
{
  switch (nucleation_parameters.model)
    {
      case Parameters::QMOMNucleationModel::constant:
        return std::make_unique<QMOMConstantNucleation>(
          nucleation_parameters.nucleation_coefficient,
          nucleation_parameters.nuclei_mean_size,
          nucleation_parameters.nuclei_standard_deviation);
      case Parameters::QMOMNucleationModel::function:
        return std::make_unique<QMOMFunctionNucleation>(
          nucleation_parameters.function_expression,
          nucleation_parameters.nuclei_mean_size,
          nucleation_parameters.nuclei_standard_deviation);
    }
  AssertThrow(false, ExcInternalError());
  return nullptr;
}

QMOMLogNormalNuclei::QMOMLogNormalNuclei(const double mean_size,
                                         const double standard_deviation)
  // log1p keeps the variance accurate for a narrow distribution
  : log_mean(std::log(mean_size) - 0.5 * std::log1p(Utilities::fixed_power<2>(
                                           standard_deviation / mean_size)))
  , log_variance(
      std::log1p(Utilities::fixed_power<2>(standard_deviation / mean_size)))
{}

void
QMOMLogNormalNuclei::fill_moment_rates(const double         nucleation_rate,
                                       std::vector<double> &moment_rates) const
{
  for (unsigned int order = 0; order < moment_rates.size(); ++order)
    moment_rates[order] =
      nucleation_rate *
      std::exp(order * log_mean + 0.5 * order * order * log_variance);
}

QMOMConstantNucleation::QMOMConstantNucleation(
  const double nucleation_coefficient,
  const double nuclei_mean_size,
  const double nuclei_standard_deviation)
  : nucleation_coefficient(nucleation_coefficient)
  , nuclei(nuclei_mean_size, nuclei_standard_deviation)
{}

void
QMOMConstantNucleation::compute_moment_rates(
  const double /*time*/,
  const std::vector<double> & /*local_scalars*/,
  std::vector<double> &moment_rates) const
{
  nuclei.fill_moment_rates(nucleation_coefficient, moment_rates);
}

QMOMFunctionNucleation::QMOMFunctionNucleation(
  const std::string &function_expression,
  const double       nuclei_mean_size,
  const double       nuclei_standard_deviation)
  : nucleation_function(1)
  , nuclei(nuclei_mean_size, nuclei_standard_deviation)
{
  // The time is the only variable: it takes the place of the spatial
  // coordinate, and the function itself does not depend on time
  nucleation_function.initialize("t",
                                 {function_expression},
                                 function_constants,
                                 false);
}

void
QMOMFunctionNucleation::compute_moment_rates(
  const double time,
  const std::vector<double> & /*local_scalars*/,
  std::vector<double> &moment_rates) const
{
  nuclei.fill_moment_rates(nucleation_function.value(Point<1>(time)),
                           moment_rates);
}

std::unique_ptr<QMOMCoalescenceModelBase>
QMOMCoalescenceModelBase::model_cast(
  const Parameters::QMOMCoalescence &coalescence_parameters)
{
  switch (coalescence_parameters.model)
    {
      case Parameters::QMOMCoalescenceModel::constant:
        return std::make_unique<QMOMConstantCoalescence>(
          coalescence_parameters.coalescence_coefficient);
      case Parameters::QMOMCoalescenceModel::function:
        return std::make_unique<QMOMFunctionCoalescence>(
          coalescence_parameters.function_expression);
    }
  AssertThrow(false, ExcInternalError());
  return nullptr;
}

void
QMOMCoalescenceModelBase::compute_moment_rates(
  const double               time,
  const std::vector<double> &local_scalars,
  const std::vector<double> &weights,
  const std::vector<double> &abscissas,
  std::vector<double>       &moment_rates) const
{
  AssertDimension(weights.size(), abscissas.size());
  std::ranges::fill(moment_rates, 0.);

  // Double sum over the ordered pairs of nodes, hence the factor 1/2. It does
  // not assume that the kernel of the function model is symmetric.
  for (unsigned int i = 0; i < abscissas.size(); ++i)
    for (unsigned int j = 0; j < abscissas.size(); ++j)
      {
        const double first_size  = abscissas[i];
        const double second_size = abscissas[j];
        if (first_size <= 0. || second_size <= 0.)
          continue;

        const double pair_rate =
          0.5 * weights[i] * weights[j] *
          compute_kernel(time, local_scalars, first_size, second_size);

        // The merged particle has the volume of the two particles
        const double merged_size =
          std::cbrt(Utilities::fixed_power<3>(first_size) +
                    Utilities::fixed_power<3>(second_size));

        // Birth of the merged particle and death of the two particles:
        // merged_size^k - first_size^k - second_size^k
        double merged_size_power = 1.;
        double first_size_power  = 1.;
        double second_size_power = 1.;
        for (double &moment_rate : moment_rates)
          {
            moment_rate += pair_rate * (merged_size_power - first_size_power -
                                        second_size_power);
            merged_size_power *= merged_size;
            first_size_power *= first_size;
            second_size_power *= second_size;
          }
      }
}

QMOMConstantCoalescence::QMOMConstantCoalescence(
  const double coalescence_coefficient)
  : coalescence_coefficient(coalescence_coefficient)
{}

double
QMOMConstantCoalescence::compute_kernel(
  const double /*time*/,
  const std::vector<double> & /*local_scalars*/,
  const double /*first_size*/,
  const double /*second_size*/) const
{
  return coalescence_coefficient;
}

QMOMFunctionCoalescence::QMOMFunctionCoalescence(
  const std::string &function_expression)
  : coalescence_function(1)
{
  // The sizes xi_1 and xi_2 of the two particles take the place of the
  // spatial coordinates
  coalescence_function.initialize("xi_1,xi_2,t",
                                  {function_expression},
                                  function_constants,
                                  true);
}

double
QMOMFunctionCoalescence::compute_kernel(
  const double time,
  const std::vector<double> & /*local_scalars*/,
  const double first_size,
  const double second_size) const
{
  coalescence_function.set_time(time);
  return coalescence_function.value(Point<2>(first_size, second_size));
}

QMOMSourceTerms::QMOMSourceTerms(const Parameters::QMOM &qmom_parameters)
{
  if (qmom_parameters.enable_growth)
    growth_model = QMOMGrowthModelBase::model_cast(qmom_parameters.growth);
  if (qmom_parameters.enable_nucleation)
    nucleation_model =
      QMOMNucleationModelBase::model_cast(qmom_parameters.nucleation);
  if (qmom_parameters.enable_coalescence)
    coalescence_model =
      QMOMCoalescenceModelBase::model_cast(qmom_parameters.coalescence);

#ifdef LETHE_WITH_QMOMKIT
  // Growth: the velocity of the nodes is their growth rate. The last
  // argument, the number of phase-space dimensions, is always 1.
  phase_space_velocity = [this](const double               time,
                                const std::vector<double> &abscissas,
                                std::vector<double>       &growth_rates,
                                const std::size_t /*n_dimensions*/) {
    if (growth_model)
      growth_model->compute_growth_rates(time,
                                         local_scalars,
                                         abscissas,
                                         growth_rates);
    else
      std::fill(growth_rates.begin(), growth_rates.end(), 0.);
  };

  // Coalescence: its rate depends on the quadrature of the current stage of
  // the integration, whose nodes are gathered from the context
  moment_source = [this](const double                      time,
                         const qmomkit::QmomSourceContext &context) {
    std::vector<double> moment_rates(context.size(), 0.);
    if (coalescence_model)
      {
        node_weights.clear();
        node_abscissas.clear();
        context.for_each_node(
          [this](const double weight, const double abscissa) {
            node_weights.push_back(weight);
            node_abscissas.push_back(abscissa);
          });
        coalescence_model->compute_moment_rates(
          time, local_scalars, node_weights, node_abscissas, moment_rates);
      }
    return qmomkit::QmomMomentRate{std::move(moment_rates)};
  };

  // Nucleation: its rate does not depend on the existing particles
  boundary_contribution = [this](const double                      time,
                                 const qmomkit::QmomSourceContext &context) {
    std::vector<double> moment_rates(context.size(), 0.);
    compute_nucleation_moment_rates(time, moment_rates);
    return qmomkit::QmomMomentRate{std::move(moment_rates)};
  };
#endif
}

void
QMOMSourceTerms::set_local_scalars(const std::vector<double> &new_local_scalars)
{
  local_scalars = new_local_scalars;
}

void
QMOMSourceTerms::compute_nucleation_moment_rates(
  const double         time,
  std::vector<double> &moment_rates) const
{
  if (nucleation_model)
    nucleation_model->compute_moment_rates(time, local_scalars, moment_rates);
  else
    std::ranges::fill(moment_rates, 0.);
}
