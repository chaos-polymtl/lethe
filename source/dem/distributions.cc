// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/distributions.h>

#include <numbers>

NormalDistribution::NormalDistribution(const double       &d_average,
                                       const double       &d_standard_deviation,
                                       const unsigned int &prn_seed,
                                       const double       &min_cutoff,
                                       const double       &max_cutoff)
  : diameter_average(d_average)
  , standard_deviation(d_standard_deviation)
  , gen(prn_seed)
{
  if (min_cutoff < 0.)
    {
      // 2.5 -> approx 99% of all diameters are bigger
      dia_min_cutoff = diameter_average - 2.5 * standard_deviation;
      AssertThrow(dia_min_cutoff > 0.,
                  ExcMessage(
                    "The \"standard deviation\" parameter is to "
                    "high relative to the \"average diameter\" "
                    "parameter. This would results with frequent negative "
                    "diameter values. To fix this problem, please change "
                    "the value of those two parameters or define a "
                    "\"minimum diameter cutoff\" bigger than \"0.\""));
    }
  else
    dia_min_cutoff = min_cutoff;

  if (max_cutoff < 0.)
    // 2.5 -> approx 99% of all diameters are smaller
    dia_max_cutoff = diameter_average + 2.5 * standard_deviation;
  else
    dia_max_cutoff = max_cutoff;

  AssertThrow(dia_min_cutoff < dia_max_cutoff,
              ExcMessage("The \"minimum diameter cutoff\" parameter need to"
                         " be smaller than the \"maximum diameter cutoff\"."));
}

void
NormalDistribution::particle_size_sampling(
  const unsigned int &number_of_particles)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(number_of_particles);

  std::normal_distribution<> distribution{diameter_average, standard_deviation};

  unsigned int n_created_diameter = 0;
  while (n_created_diameter < number_of_particles)
    {
      const double temp_diameter = distribution(gen);
      if (temp_diameter > dia_min_cutoff && temp_diameter < dia_max_cutoff)
        {
          n_created_diameter++;
          this->particle_sizes.push_back(temp_diameter);
        }
    }
}

double
NormalDistribution::find_min_diameter()
{
  return dia_min_cutoff;
}

double
NormalDistribution::find_max_diameter()
{
  return dia_max_cutoff;
}

void
NormalDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  const double z_min = (dia_min_cutoff - diameter_average) / standard_deviation;
  const double z_max = (dia_max_cutoff - diameter_average) / standard_deviation;

  // Quantile
  const double q_min = 0.5 * (1.0 + std::erf(z_min / std::numbers::sqrt2));
  const double q_max = 0.5 * (1.0 + std::erf(z_max / std::numbers::sqrt2));

  pcout << "The particle size distribution of particle type " << particle_type
        << " is normal." << std::endl;

  if (q_min > 0.25 || q_max < 0.75)
    pcout
      << "Warning: The minimal and maximal cutoffs of the distribution are at "
      << q_min * 100. << " and " << q_max * 100. << " respectively."
      << std::endl;
  else
    pcout << "The minimal and maximal cutoffs of the distribution are at "
          << q_min * 100. << " and " << q_max * 100. << " respectively."
          << std::endl;
}

LogNormalDistribution::LogNormalDistribution(const double &d_average,
                                             const double &d_standard_deviation,
                                             const unsigned int &prn_seed,
                                             const double        min_cutoff,
                                             const double        max_cutoff)
  : sigma_ln(std::sqrt(std::log(
      1. + Utilities::fixed_power<2>(d_standard_deviation / d_average))))
  , mu_ln(std::log(d_average) -
          0.5 * Utilities::fixed_power<2>(d_standard_deviation))
  , gen(prn_seed)
{
  if (min_cutoff < 0.)
    {
      // approx 99% of all diameters are bigger
      dia_min_cutoff = std::exp(mu_ln - 2.5 * sigma_ln);
      AssertThrow(dia_min_cutoff > 0.,
                  ExcMessage(
                    "The \"standard deviation\" parameter is to "
                    "high relative to the \"average diameter\" "
                    "parameter. This would results with frequent negative "
                    "diameter values. To fix this problem, please change "
                    "the value of those two parameters or define a "
                    "\"minimum diameter cutoff\" bigger than \"0.\""));
    }
  else
    dia_min_cutoff = min_cutoff;

  if (max_cutoff < 0.)
    // approx 99% of all diameters are bigger
    dia_max_cutoff = std::exp(mu_ln + 2.5 * sigma_ln);
  else
    dia_max_cutoff = max_cutoff;

  AssertThrow(dia_min_cutoff < dia_max_cutoff,
              ExcMessage(
                "The \"minimum diameter cutoff\" parameter need to be smaller "
                "than the \"maximum diameter cutoff\"."));
}


void
LogNormalDistribution::particle_size_sampling(
  const unsigned int &number_of_particles)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(number_of_particles);

  std::lognormal_distribution<> distribution{mu_ln, sigma_ln};

  unsigned int n_created_diameter = 0;
  while (n_created_diameter < number_of_particles)
    {
      const double temp_diameter = distribution(gen);
      if (temp_diameter > dia_min_cutoff && temp_diameter < dia_max_cutoff)
        {
          n_created_diameter++;
          this->particle_sizes.push_back(temp_diameter);
        }
    }
}

double
LogNormalDistribution::find_min_diameter()
{
  return dia_min_cutoff;
}

double
LogNormalDistribution::find_max_diameter()
{
  return dia_max_cutoff;
}

void
LogNormalDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  const double z_min = (std::log(dia_min_cutoff) - mu_ln) / sigma_ln;
  const double z_max = (std::log(dia_max_cutoff) - mu_ln) / sigma_ln;

  // Quantile
  const double q_min = 0.5 * (1.0 + std::erf(z_min / std::numbers::sqrt2));
  const double q_max = 0.5 * (1.0 + std::erf(z_max / std::numbers::sqrt2));

  pcout << "The particle size distribution of particle type " << particle_type
        << " is lognormal." << std::endl;

  if (q_min > 0.25 || q_max < 0.75)
    pcout
      << "Warning: The minimal and maximal cutoffs of the distribution are at "
      << q_min * 100. << " and " << q_max * 100. << " percentiles respectively."
      << std::endl;
  else
    pcout << "The minimal and maximal cutoffs of the distribution are at "
          << q_min * 100. << " and " << q_max * 100.
          << " percentiles respectively." << std::endl;
}

UniformDistribution::UniformDistribution(const double &d_values)
  : diameter_value(d_values)
{}

void
UniformDistribution::particle_size_sampling(
  const unsigned int &number_of_particles)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(number_of_particles);

  for (unsigned int n = 0; n < number_of_particles; ++n)
    this->particle_sizes.push_back(this->diameter_value);
}

double
UniformDistribution::find_min_diameter()
{
  return this->diameter_value;
}

double
UniformDistribution::find_max_diameter()
{
  return this->diameter_value;
}

void
UniformDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  pcout << "The particle size distribution of particle type " << particle_type
        << " is uniform." << std::endl;
}

CustomDistribution::CustomDistribution(
  const std::vector<double> &d_list,
  const std::vector<double> &d_probabilities,
  const unsigned int        &prn_seed)
  : diameter_custom_values(d_list)
  , gen(prn_seed)
{
  std::vector<double> cumulative_probability_vector, n_i_vector;
  double              n_tot = 0.;

  cumulative_probability_vector.reserve(d_probabilities.size());
  n_i_vector.reserve(d_probabilities.size());

  for (unsigned int i = 0; i < diameter_custom_values.size(); ++i)
    {
      n_i_vector.push_back(d_probabilities[i] / Utilities::fixed_power<3>(
                                                  diameter_custom_values[i]));
      n_tot += n_i_vector[i];
    }

  double i_probability;
  double cumulative_value = 0.;

  for (unsigned int i = 0; i < diameter_custom_values.size(); ++i)
    {
      i_probability = n_i_vector[i] / n_tot;
      cumulative_value += i_probability;
      cumulative_probability_vector.push_back(cumulative_value);
    }

  diameter_custom_cumul_prob = cumulative_probability_vector;
}

void
CustomDistribution::particle_size_sampling(
  const unsigned int &number_of_particles)
{
  this->particle_sizes.clear();
  this->particle_sizes.reserve(number_of_particles);

  std::uniform_real_distribution<> dis(0.0, diameter_custom_cumul_prob.back());

  for (unsigned int i = 0; i < number_of_particles; ++i)
    {
      // Search to find the appropriate diameter index
      auto it = std::ranges::upper_bound(diameter_custom_cumul_prob, dis(gen));

      // if dis(gen) returns exactly the maximum value of the cumulative
      // distribution vector
      if (it == diameter_custom_cumul_prob.end())
        {
          it = it - 1; // Move back to the last valid element
        }

      unsigned int index =
        static_cast<unsigned int>(it - diameter_custom_cumul_prob.begin());

      this->particle_sizes.push_back(diameter_custom_values[index]);
    }
}

double
CustomDistribution::find_min_diameter()
{
  return *std::ranges::min_element(diameter_custom_values);
}

double
CustomDistribution::find_max_diameter()
{
  return *std::ranges::max_element(diameter_custom_values);
}

void
CustomDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  pcout << "The particle size distribution of particle type " << particle_type
        << " is custom." << std::endl;
}
