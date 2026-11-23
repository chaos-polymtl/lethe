// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/parameters_multiphysics.h>
#include <core/utilities.h>

#include <solvers/qmom_source_terms.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <magic_enum.hpp>

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
  AssertThrow(false,
              ExcMessage("The nucleation model \"" +
                         enum_to_string(nucleation_parameters.model) +
                         "\" of the QMOM physics is not handled by "
                         "QMOMNucleationModelBase::model_cast(). Add a case "
                         "for this model to its switch."));
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
  AssertThrow(false,
              ExcMessage("The coalescence model \"" +
                         enum_to_string(coalescence_parameters.model) +
                         "\" of the QMOM physics is not handled by "
                         "QMOMCoalescenceModelBase::model_cast(). Add a case "
                         "for this model to its switch."));
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

QMOMLipidNanoparticleMixture::QMOMLipidNanoparticleMixture(
  const Parameters::QMOMLipidNanoparticleProduction &process_parameters,
  const Parameters::QMOMPhysicalConstants           &physical_constants)
  : parameters(process_parameters)
  , constants(physical_constants)
  , thermal_energy(physical_constants.boltzmann_constant *
                   process_parameters.operating_temperature)
  , capillary_length(4. * process_parameters.surface_tension *
                     process_parameters.lipid_molecular_volume / thermal_energy)
{}

QMOMLipidNanoparticleLocalState
QMOMLipidNanoparticleMixture::compute_local_state(
  const std::vector<double> &local_scalars) const
{
  AssertDimension(
    local_scalars.size(),
    magic_enum::enum_count<Parameters::QMOMLipidNanoparticleScalar>());

  // The transported concentrations may slightly overshoot their physical
  // bounds, so the volume fraction is clipped to [0, 1] and the free lipid
  // concentration to non-negative values
  const double organic_volume_fraction = std::clamp(
    local_scalars[Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::ethanol_concentration)] /
      parameters.density_organic,
    0.,
    1.);
  const double aqueous_volume_fraction = 1. - organic_volume_fraction;

  // ln x = phi_a ln x_a + phi_o ln x_o
  const double mole_fraction_solubility =
    std::exp(aqueous_volume_fraction *
               std::log(parameters.mole_fraction_solubility_aqueous) +
             organic_volume_fraction *
               std::log(parameters.mole_fraction_solubility_organic));

  // C_inf = M_l x / (1 - x) (phi_a rho_a / M_a + phi_o rho_o / M_o)
  const double solubility =
    parameters.lipid_molar_mass * mole_fraction_solubility /
    (1. - mole_fraction_solubility) *
    (aqueous_volume_fraction * parameters.density_aqueous /
       parameters.molar_mass_aqueous +
     organic_volume_fraction * parameters.density_organic /
       parameters.molar_mass_organic);

  const double free_lipids_concentration = std::max(
    local_scalars[Parameters::get_scalar_index(
      Parameters::QMOMLipidNanoparticleScalar::free_lipids_concentration)],
    0.);

  return {aqueous_volume_fraction,
          solubility,
          free_lipids_concentration / solubility};
}

double
QMOMLipidNanoparticleMixture::compute_equilibrium_supersaturation(
  const double size) const
{
  Assert(size > 0., ExcLowerRangeType<double>(size, 0.));
  return std::exp(capillary_length / size);
}

double
QMOMLipidNanoparticleMixture::compute_critical_size(
  const double supersaturation) const
{
  Assert(supersaturation > 1., ExcLowerRangeType<double>(supersaturation, 1.));
  return capillary_length / std::log(supersaturation);
}

double
QMOMLipidNanoparticleMixture::compute_inverse_debye_length(
  const double aqueous_volume_fraction) const
{
  // The permittivity depends on the local composition of the mixture
  const double relative_permittivity =
    aqueous_volume_fraction * parameters.permittivity_aqueous +
    (1. - aqueous_volume_fraction) * parameters.permittivity_organic;

  return std::sqrt(
    2. * Utilities::fixed_power<2>(constants.elementary_charge) *
    constants.avogadro_number * parameters.ion_concentration_mixture /
    (constants.permittivity_vacuum * relative_permittivity * thermal_energy));
}

double
QMOMLipidNanoparticleMixture::compute_electrostatic_prefactor(
  const double aqueous_volume_fraction) const
{
  const double relative_permittivity =
    aqueous_volume_fraction * parameters.permittivity_aqueous +
    (1. - aqueous_volume_fraction) * parameters.permittivity_organic;
  const double potential_factor =
    std::tanh(constants.elementary_charge * parameters.surface_potential /
              (4. * thermal_energy));

  // 128 pi N_A c_i k_B T / kappa^2, with kappa^2 proportional to c_i, is
  // written without c_i, so that it stays finite without ions
  return 64. * numbers::PI * constants.permittivity_vacuum *
         relative_permittivity *
         Utilities::fixed_power<2>(thermal_energy /
                                   constants.elementary_charge) *
         Utilities::fixed_power<2>(potential_factor);
}

QMOMLipidNanoparticleGrowth::QMOMLipidNanoparticleGrowth(
  const QMOMLipidNanoparticleMixture &mixture)
  : mixture(mixture)
{}

void
QMOMLipidNanoparticleGrowth::compute_growth_rates(
  const double /*time*/,
  const std::vector<double> &local_scalars,
  const std::vector<double> &abscissas,
  std::vector<double>       &growth_rates) const
{
  AssertDimension(growth_rates.size(), abscissas.size());

  const QMOMLipidNanoparticleLocalState state =
    mixture.compute_local_state(local_scalars);

  for (unsigned int node = 0; node < abscissas.size(); ++node)
    {
      const double size = abscissas[node];
      if (size <= 0.)
        {
          growth_rates[node] = 0.;
          continue;
        }

      // G = k_g xi C_inf (S - S*(xi)), negative below the equilibrium size
      growth_rates[node] = mixture.parameters.growth_factor * size *
                           state.solubility *
                           (state.supersaturation -
                            mixture.compute_equilibrium_supersaturation(size));
    }
}

QMOMLipidNanoparticleNucleation::QMOMLipidNanoparticleNucleation(
  const QMOMLipidNanoparticleMixture &mixture)
  : mixture(mixture)
  , exponent_coefficient(
      16. * numbers::PI *
      Utilities::fixed_power<3>(mixture.parameters.surface_tension) *
      Utilities::fixed_power<2>(mixture.parameters.lipid_molecular_volume) /
      (3. * Utilities::fixed_power<3>(mixture.thermal_energy)))
{}

double
QMOMLipidNanoparticleNucleation::compute_nucleation_rate(
  const double supersaturation) const
{
  // Nucleation only happens in a supersaturated mixture
  if (supersaturation <= 1.)
    return 0.;

  return mixture.parameters.nucleation_factor *
         std::exp(-exponent_coefficient /
                  Utilities::fixed_power<2>(std::log(supersaturation)));
}

void
QMOMLipidNanoparticleNucleation::compute_moment_rates(
  const double /*time*/,
  const std::vector<double> &local_scalars,
  std::vector<double>       &moment_rates) const
{
  const double supersaturation =
    mixture.compute_local_state(local_scalars).supersaturation;
  const double nucleation_rate = compute_nucleation_rate(supersaturation);

  // Close to saturation, the critical size diverges while the rate vanishes:
  // the moments of the nuclei are not evaluated, to avoid 0 * inf
  if (!(nucleation_rate > 0.))
    {
      std::ranges::fill(moment_rates, 0.);
      return;
    }

  const QMOMLogNormalNuclei nuclei(
    mixture.compute_critical_size(supersaturation),
    mixture.parameters.nuclei_standard_deviation);
  nuclei.fill_moment_rates(nucleation_rate, moment_rates);
}

QMOMLipidNanoparticleCoalescence::QMOMLipidNanoparticleCoalescence(
  const QMOMLipidNanoparticleMixture &mixture)
  : mixture(mixture)
  , log_step(std::log(max_separation / min_separation) / n_intervals)
  , separations(n_intervals + 1)
{
  static_assert(n_intervals % 2 == 0,
                "The Simpson rule requires an even number of intervals.");

  // The end points are set exactly, so that the integral covers the whole
  // interval despite the round-off of the exponential
  for (unsigned int i = 0; i <= n_intervals; ++i)
    separations[i] = min_separation * std::exp(i * log_step);
  separations.front() = min_separation;
  separations.back()  = max_separation;
}

double
QMOMLipidNanoparticleCoalescence::compute_stability_ratio(
  const double first_size,
  const double second_size,
  const double aqueous_volume_fraction) const
{
  const double first_radius   = 0.5 * first_size;
  const double second_radius  = 0.5 * second_size;
  const double radius_sum     = first_radius + second_radius;
  const double radius_product = first_radius * second_radius;

  // Terms of the electrostatic repulsion that do not depend on the distance:
  // Phi_el(a) = electrostatic_scale exp(-kappa a)
  const double inverse_debye_length =
    mixture.compute_inverse_debye_length(aqueous_volume_fraction);
  const double electrostatic_scale =
    mixture.compute_electrostatic_prefactor(aqueous_volume_fraction) *
    radius_product / radius_sum;

  const double hamaker_constant = mixture.parameters.hamaker_constant;
  const double thermal_energy   = mixture.thermal_energy;

  // Integrand in the logarithm of the distance a, hence the factor a of
  // da = a d(ln a): exp(Phi(a) / (k_B T)) a / (a + r_s)^2
  const auto integrand = [&](const double separation) {
    // van der Waals attraction of two spheres
    const double near_term =
      Utilities::fixed_power<2>(separation) + 2. * separation * radius_sum;
    const double far_term = near_term + 4. * radius_product;
    const double van_der_waals_energy =
      -hamaker_constant / 6. *
      (2. * radius_product / near_term + 2. * radius_product / far_term +
       std::log(near_term / far_term));

    const double electrostatic_energy =
      electrostatic_scale * std::exp(-inverse_debye_length * separation);

    return std::exp((van_der_waals_energy + electrostatic_energy) /
                    thermal_energy) *
           separation / Utilities::fixed_power<2>(separation + radius_sum);
  };

  // Simpson rule on the logarithmic grid
  double integral =
    integrand(separations.front()) + integrand(separations.back());
  for (unsigned int i = 1; i < n_intervals; ++i)
    integral += (i % 2 == 0 ? 2. : 4.) * integrand(separations[i]);
  integral *= log_step / 3.;

  // Tail beyond max_separation, where the interaction energy vanishes
  integral += 1. / (max_separation + radius_sum);

  return radius_sum * integral;
}

double
QMOMLipidNanoparticleCoalescence::compute_kernel(
  const double /*time*/,
  const std::vector<double> &local_scalars,
  const double               first_size,
  const double               second_size) const
{
  const Parameters::QMOMLipidNanoparticleProduction &parameters =
    mixture.parameters;
  const double aqueous_volume_fraction =
    mixture.compute_local_state(local_scalars).aqueous_volume_fraction;

  // Brownian (Smoluchowski) kernel
  const double brownian_kernel =
    2. * mixture.thermal_energy *
    Utilities::fixed_power<2>(first_size + second_size) /
    (3. * parameters.viscosity_mixture * first_size * second_size);

  // Attachment efficiency of the Fuchs model, limited by the PEG lipids
  const double attachment_efficiency = std::min(
    parameters.peg_factor /
      compute_stability_ratio(first_size, second_size, aqueous_volume_fraction),
    parameters.peg_factor);

  return brownian_kernel * attachment_efficiency;
}

QMOMSourceTerms::QMOMSourceTerms(const Parameters::QMOM &qmom_parameters)
{
  switch (qmom_parameters.process)
    {
      case Parameters::QMOMProcess::generic:
        if (qmom_parameters.enable_growth)
          growth_model =
            QMOMGrowthModelBase::model_cast(qmom_parameters.growth);
        if (qmom_parameters.enable_nucleation)
          nucleation_model =
            QMOMNucleationModelBase::model_cast(qmom_parameters.nucleation);
        if (qmom_parameters.enable_coalescence)
          coalescence_model =
            QMOMCoalescenceModelBase::model_cast(qmom_parameters.coalescence);
        break;
      case Parameters::QMOMProcess::lipid_nanoparticle_production:
        // The models of the process replace those of the "growth",
        // "nucleation" and "coalescence" subsections, but are still enabled
        // one by one
        lipid_mixture = std::make_unique<QMOMLipidNanoparticleMixture>(
          qmom_parameters.lipid_nanoparticle_production,
          qmom_parameters.physical_constants);
        if (qmom_parameters.enable_growth)
          growth_model =
            std::make_unique<QMOMLipidNanoparticleGrowth>(*lipid_mixture);
        if (qmom_parameters.enable_nucleation)
          nucleation_model =
            std::make_unique<QMOMLipidNanoparticleNucleation>(*lipid_mixture);
        if (qmom_parameters.enable_coalescence)
          coalescence_model =
            std::make_unique<QMOMLipidNanoparticleCoalescence>(*lipid_mixture);
        break;
    }

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
      std::ranges::fill(growth_rates, 0.);
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

bool
QMOMSourceTerms::update_local_scalars(const std::vector<double> &moments_before,
                                      const std::vector<double> &moments_after,
                                      std::vector<double>       &scalars) const
{
  if (!lipid_mixture)
    return true;

  // Order of the moment proportional to the volume of the particles
  constexpr unsigned int volume_moment_order = 3;
  AssertIndexRange(volume_moment_order, moments_before.size());
  AssertDimension(moments_before.size(), moments_after.size());
  AssertDimension(
    scalars.size(),
    magic_enum::enum_count<Parameters::QMOMLipidNanoparticleScalar>());

  // The mass of lipid per unit volume in the particles is rho_l k_v M_3
  const Parameters::QMOMLipidNanoparticleProduction &parameters =
    lipid_mixture->parameters;
  const double lipid_mass_per_volume_moment =
    parameters.lipid_density * parameters.shape_factor;
  double      &free_lipids_concentration = scalars[Parameters::get_scalar_index(
    Parameters::QMOMLipidNanoparticleScalar::free_lipids_concentration)];
  const double free_lipids_concentration_before = free_lipids_concentration;
  free_lipids_concentration -=
    lipid_mass_per_volume_moment *
    (moments_after[volume_moment_order] - moments_before[volume_moment_order]);

  if (free_lipids_concentration >= 0.)
    return true;

  // Round-off of the lipid mass, free and in the particles
  constexpr double round_off_relative_tolerance = 1e-12;
  const double     lipid_mass_scale =
    std::abs(free_lipids_concentration_before) +
    lipid_mass_per_volume_moment *
      std::max(std::abs(moments_before[volume_moment_order]),
               std::abs(moments_after[volume_moment_order]));
  if (free_lipids_concentration >=
      -round_off_relative_tolerance * lipid_mass_scale)
    {
      free_lipids_concentration = 0.;
      return true;
    }
  return false;
}
