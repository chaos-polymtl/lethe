// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/qmom_field_function.h>

#include <deal.II/base/exceptions.h>

#include <string>
#include <vector>

#ifdef LETHE_WITH_QMOMKIT
#  include <qmomkit/distributions.hpp>
#endif

template <int dim>
QMOMFieldFunction<dim>::QMOMFieldFunction(
  const Parameters::QMOMFunctions<dim> &functions,
  const unsigned int                    n_scalars,
  const unsigned int                    n_moments,
  const std::string                    &description)
  : Function<dim>(n_scalars + n_moments)
  , functions(functions)
  , n_scalars(n_scalars)
  , n_moments(n_moments)
  , scalars_are_broadcast(functions.scalars->n_components == 1)
{
  AssertThrow(
    scalars_are_broadcast || functions.scalars->n_components == n_scalars,
    ExcMessage("The scalars function of the QMOM " + description + " has " +
               std::to_string(functions.scalars->n_components) +
               " expressions, but the QMOM physics transports " +
               std::to_string(n_scalars) +
               " scalars (the \"number of scalars\" of the QMOM "
               "subsection). Give either a single expression, applied to "
               "every scalar, or one expression per scalar, separated by "
               "semicolons."));

#ifndef LETHE_WITH_QMOMKIT
  // The moments of a particle size distribution are computed by QMOMKit
  AssertThrow(functions.distribution_type ==
                Parameters::QMOMDistributionType::none,
              QMOMWithoutQMOMKitError());
#endif
}

template <int dim>
double
QMOMFieldFunction<dim>::value(const Point<dim>  &p,
                              const unsigned int component) const
{
  AssertIndexRange(component, this->n_components);

  if (component < n_scalars)
    return functions.scalars->value(p, scalars_are_broadcast ? 0 : component);

  return compute_moments(p)[component - n_scalars];
}

template <int dim>
void
QMOMFieldFunction<dim>::vector_value(const Point<dim> &p,
                                     Vector<double>   &values) const
{
  AssertDimension(values.size(), this->n_components);

  for (unsigned int i_scalar = 0; i_scalar < n_scalars; ++i_scalar)
    values[i_scalar] =
      functions.scalars->value(p, scalars_are_broadcast ? 0 : i_scalar);

  const std::vector<double> moments = compute_moments(p);
  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    values[n_scalars + i_moment] = moments[i_moment];
}

template <int dim>
void
QMOMFieldFunction<dim>::set_time(const double new_time)
{
  Function<dim>::set_time(new_time);
  functions.set_time(new_time);
}

template <int dim>
std::vector<double>
QMOMFieldFunction<dim>::compute_moments(
  [[maybe_unused]] const Point<dim> &p) const
{
  using DistributionType = Parameters::QMOMDistributionType;

  std::vector<double> moments(n_moments, 0.);
  if (functions.distribution_type == DistributionType::none)
    return moments;

#ifdef LETHE_WITH_QMOMKIT
  // Number density followed by the parameters of the distribution
  Vector<double> distribution_values(functions.distribution->n_components);
  functions.distribution->vector_value(p, distribution_values);
  const double number_density = distribution_values[0];

  AssertThrow(number_density >= 0.,
              ExcMessage("The number density of the particle size "
                         "distribution of the QMOM physics must be positive "
                         "or zero, but it is " +
                         std::to_string(number_density) + " particles/m^3."));

  // Without particles, all the moments are zero whatever the distribution
  if (number_density <= 0.)
    return moments;

  qmomkit::DistributionSpecification distribution;
  switch (functions.distribution_type)
    {
      case DistributionType::normal:
        distribution =
          qmomkit::NormalDistribution{.mean = distribution_values[1],
                                      .standard_deviation =
                                        distribution_values[2]};
        break;
      case DistributionType::lognormal:
        distribution =
          qmomkit::LogNormalDistribution{.normal_mean = distribution_values[1],
                                         .normal_standard_deviation =
                                           distribution_values[2]};
        break;
      case DistributionType::gamma:
        distribution =
          qmomkit::GammaDistribution{.shape = distribution_values[1],
                                     .scale = distribution_values[2]};
        break;
      case DistributionType::weibull:
        distribution =
          qmomkit::WeibullDistribution{.shape = distribution_values[1],
                                       .scale = distribution_values[2]};
        break;
      case DistributionType::beta:
        distribution =
          qmomkit::BetaDistribution{.alpha = distribution_values[1],
                                    .beta  = distribution_values[2]};
        break;
      case DistributionType::chi:
        distribution = qmomkit::ChiDistribution{.degrees_of_freedom =
                                                  distribution_values[1]};
        break;
      case DistributionType::chi_squared:
        distribution = qmomkit::ChiSquaredDistribution{
          .degrees_of_freedom = distribution_values[1]};
        break;
      case DistributionType::none:
        Assert(false, ExcInternalError());
        break;
    }

  const qmomkit::RawMomentResult result =
    qmomkit::compute_raw_moments(distribution, n_moments - 1, number_density);
  AssertThrow(result.ok(),
              ExcMessage("QMOMKit could not compute the moments of the "
                         "particle size distribution of the QMOM physics: " +
                         result.message +
                         " Check the parameters of the distribution."));

  for (unsigned int i_moment = 0; i_moment < n_moments; ++i_moment)
    moments[i_moment] = result.moments[i_moment];
#endif

  return moments;
}

template class QMOMFieldFunction<2>;
template class QMOMFieldFunction<3>;
