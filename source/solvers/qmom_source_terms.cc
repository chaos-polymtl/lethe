// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_multiphysics.h>

#include <solvers/qmom_source_terms.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

QMOMPhaseSpaceVelocity
no_growth_velocity()
{
  return [](const double /*time*/,
            const std::vector<double> & /*abscissas*/,
            std::vector<double> &growth_rates,
            const std::size_t /*n_phase_space_dimensions*/) {
    std::fill(growth_rates.begin(), growth_rates.end(), 0.);
  };
}

QMOMPhaseSpaceVelocity
constant_growth_velocity(const double growth_coefficient)
{
  return [growth_coefficient](const double /*time*/,
                              const std::vector<double> & /*abscissas*/,
                              std::vector<double> &growth_rates,
                              const std::size_t /*n_phase_space_dimensions*/) {
    std::fill(growth_rates.begin(), growth_rates.end(), growth_coefficient);
  };
}

QMOMPhaseSpaceVelocity
linear_growth_velocity(const double growth_coefficient)
{
  return [growth_coefficient](const double /*time*/,
                              const std::vector<double> &abscissas,
                              std::vector<double>       &growth_rates,
                              const std::size_t /*n_phase_space_dimensions*/) {
    AssertDimension(growth_rates.size(), abscissas.size());
    for (unsigned int node = 0; node < abscissas.size(); ++node)
      growth_rates[node] = growth_coefficient * abscissas[node];
  };
}

QMOMPhaseSpaceVelocity
function_growth_velocity(const std::string &function_expression)
{
  // The variable of the function is the particle size xi instead of a
  // spatial coordinate. The parser is shared by the copies of the returned
  // velocity.
  auto growth_function = std::make_shared<FunctionParser<1>>(1);
  const std::map<std::string, double> constants{{"pi", numbers::PI}};
  growth_function->initialize("xi,t", {function_expression}, constants, true);

  return [growth_function](const double               time,
                           const std::vector<double> &abscissas,
                           std::vector<double>       &growth_rates,
                           const std::size_t /*n_phase_space_dimensions*/) {
    AssertDimension(growth_rates.size(), abscissas.size());
    growth_function->set_time(time);
    for (unsigned int node = 0; node < abscissas.size(); ++node)
      growth_rates[node] = growth_function->value(Point<1>(abscissas[node]));
  };
}

QMOMSourceTerms::QMOMSourceTerms(const Parameters::QMOM &qmom_parameters)
  : enable_growth(qmom_parameters.enable_growth)
{
  if (!enable_growth)
    {
      growth_velocity = no_growth_velocity();
      return;
    }

  const Parameters::QMOMGrowth &growth = qmom_parameters.growth;
  switch (growth.model)
    {
      case Parameters::QMOMGrowthModel::constant:
        growth_velocity = constant_growth_velocity(growth.growth_coefficient);
        break;
      case Parameters::QMOMGrowthModel::linear:
        growth_velocity = linear_growth_velocity(growth.growth_coefficient);
        break;
      case Parameters::QMOMGrowthModel::function:
        growth_velocity = function_growth_velocity(growth.function_expression);
        break;
    }
}
