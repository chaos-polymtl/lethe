// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/distributions.h>

#include <deal.II/lac/lapack_full_matrix.h>

#include <numbers>
Distribution::Distribution(
  const DistributionWeightingType &distribution_weighting_type)
  : weighting_type(distribution_weighting_type)
{}

NormalDistribution::NormalDistribution(
  const double                    &d_average,
  const double                    &d_standard_deviation,
  const unsigned int              &prn_seed,
  const double                     min_cutoff,
  const double                     max_cutoff,
  const DistributionWeightingType &distribution_weighting_type)
  : Distribution(distribution_weighting_type)
  , gen(prn_seed)
{
  // We need the average and the standard deviation to be based on number to do
  // our sampling.
  if (this->weighting_type == DistributionWeightingType::number_based)
    {
      diameter_average   = d_average;
      standard_deviation = d_standard_deviation;
    }
  // If the user specified the distribution as volume-based, we need to convert
  // the average and standard deviation to be number-based. To do so, we need to
  // solve a nonlinear system.
  // https://mfix.netl.doe.gov/doc/mfix-exa/guide/latest/references/size_distributions.html#normal-distribution
  else if (this->weighting_type == DistributionWeightingType::volume_based)
    {
      // Either there is an error with the Newton-Raphson solver or the problem
      // is too stiff due to the order of the polynomial we are trying to solve
      // for. Because of this, volume-based distribution is not supported for
      // now. MFIX-exa uses a homotopy method to solve the system.

      // Volume-based parameters
      const double mv  = d_average; // volume-based mean
      const double mv2 = mv * mv;   // volume-based mean squared

      const double sv = d_standard_deviation; // volume-based standard deviation
      const double sv2 = sv * sv; // volume-based standard deviation squared

      // Matrix and RHS
      LAPACKFullMatrix<double> system_matrix;
      Vector<double>           system_rhs;

      double       residual_l2_norm  = 1.;
      unsigned int iteration_counter = 0;

      // Initial estimates
      // mean based on number
      double mu_n = d_average;
      // standard deviation based on number
      double sigma_n = d_standard_deviation;

      system_matrix.reinit(2);
      system_rhs.reinit(2);
      while (residual_l2_norm > 1.e-12 && iteration_counter < 100)
        {
          const double mn  = mu_n;      // mean
          const double mn2 = mn * mn;   // mean squared
          const double mn3 = mn2 * mn;  // mean cubed
          const double mn4 = mn2 * mn2; // mean ^4
          const double mn5 = mn4 * mn;  // mean ^5

          const double sn  = sigma_n;   // standard deviation
          const double sn2 = sn * sn;   // standard deviation squared
          const double sn3 = sn2 * sn;  // standard deviation cubed
          const double sn4 = sn2 * sn2; // standard deviation ^4

          // Jacobian
          const double dR1_sigma_n =
            6. * (2. * (sn3 + sn * mn2) - mv * sn * mn);

          const double dR1_mu_n =
            12. * sn2 * mn + 4 * mn3 - 3. * mv * (sn2 + mn2);

          const double dR2_sigma_n =
            12. * sn3 * (5. * mn - 2. * mv) +
            2. * sn * (10. * mn3 - 12. * mn2 * mv + 3. * mn * mv2) -
            6. * sv2 * sn * mn;

          const double dR2_mu_n = 15. * sn4 +
                                  sn2 * (30. * mn2 - 24. * mn * mv + 3. * mv2) +
                                  5. * mn4 + 3. * mn2 * mv2 - 8. * mn3 * mv -
                                  sv2 * (3. * sn2 + 3. * mn2);

          system_matrix.set(0, 0, dR1_sigma_n);
          system_matrix.set(0, 1, dR1_mu_n);
          system_matrix.set(1, 0, dR2_sigma_n);
          system_matrix.set(1, 1, dR2_mu_n);

          // Resisual RHS
          const double R1 =
            3. * sn4 + 6. * sn2 * mn2 + mn4 - mv * (3. * sn2 * mn + mn3);

          const double R2 = 3. * sn4 * (5. * mn - 2. * mv) +
                            sn2 * (10. * mn3 - 12. * mn2 * mv + 3. * mn * mv2) +
                            mn5 + mn3 * mv2 - 2. * mn4 * mv -
                            sv2 * (3. * sn2 * mn + mn3);

          // Add the minus sign for Newton-Raphson
          system_rhs[0] = -R1;
          system_rhs[1] = -R2;

          // Solve for the update
          system_matrix.compute_lu_factorization();
          system_matrix.solve(system_rhs);

          // USE THE MATHEMATICALLY CORRECT MAPPING:
          sigma_n += system_rhs[0];
          mu_n += system_rhs[1];

          // While loop conditions
          residual_l2_norm = system_rhs.l2_norm();
          iteration_counter++;
        }
      AssertThrow(
        iteration_counter < 100.,
        ExcMessage("The numbered weighted mean and standard deviation of "
                   "the normal distribution has not been found from the volume "
                   "weighted values provided in the parameter file. To solve "
                   "this issue, define the distribution as number-based."));

      diameter_average   = mu_n;
      standard_deviation = sigma_n;
    }

  if (min_cutoff < 0.)
    {
      // 2.5 -> approx 99% of all diameters are bigger
      this->dia_min_cutoff = diameter_average - 2.5 * standard_deviation;
      AssertThrow(this->dia_min_cutoff > 0.,
                  ExcMessage(
                    "The \"standard deviation\" parameter is too "
                    "high relative to the \"average diameter\" "
                    "parameter. This would results with frequent negative "
                    "diameter values. To fix this problem, please change "
                    "the value of those two parameters or define a "
                    "\"minimum diameter cutoff\" bigger than \"0.\""));
    }
  else
    this->dia_min_cutoff = min_cutoff;

  if (max_cutoff < 0.)
    // 2.5 -> approx 99% of all diameters are smaller
    this->dia_max_cutoff = diameter_average + 2.5 * standard_deviation;
  else
    this->dia_max_cutoff = max_cutoff;

  AssertThrow(this->dia_min_cutoff < this->dia_max_cutoff,
              ExcMessage("The \"minimum diameter cutoff\" parameter needs to"
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
      if (temp_diameter > this->dia_min_cutoff &&
          temp_diameter < this->dia_max_cutoff)
        {
          n_created_diameter++;
          this->particle_sizes.push_back(temp_diameter);
        }
    }
}

double
NormalDistribution::find_min_diameter()
{
  return this->dia_min_cutoff;
}

double
NormalDistribution::find_max_diameter()
{
  return this->dia_max_cutoff;
}

void
NormalDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  const double z_min =
    (this->dia_min_cutoff - diameter_average) / standard_deviation;
  const double z_max =
    (this->dia_max_cutoff - diameter_average) / standard_deviation;

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

LogNormalDistribution::LogNormalDistribution(
  const double       &d_average,
  const double       &d_standard_deviation,
  const unsigned int &prn_seed,
  double              min_cutoff,
  double              max_cutoff,
  const Parameters::Lagrangian::DistributionWeightingType
    &distribution_weighting_type)
  : Distribution(distribution_weighting_type)
  , sigma_ln(std::sqrt(std::log(
      1. + Utilities::fixed_power<2>(d_standard_deviation / d_average))))
  , gen(prn_seed)
{
  // We need the number-based average and standard deviation to do
  // our sample. The standard deviation is always the same.
  if (this->weighting_type == DistributionWeightingType::number_based)
    {
      mu_ln = std::log(d_average) - 0.5 * Utilities::fixed_power<2>(sigma_ln);
    }
  // If the user specified the distribution as volume-based, we need to convert
  // the average and standard deviation to be number-based.
  // https://mfix.netl.doe.gov/doc/mfix-exa/guide/latest/references/size_distributions.html#log-normal-distribution
  else if (this->weighting_type == DistributionWeightingType::volume_based)
    {
      mu_ln = std::log(d_average) - 3.5 * Utilities::fixed_power<2>(sigma_ln);
    }

  if (min_cutoff < 0.)
    {
      // approx 99% of all diameters are bigger
      this->dia_min_cutoff = std::exp(mu_ln - 2.5 * sigma_ln);
      AssertThrow(this->dia_min_cutoff > 0.,
                  ExcMessage(
                    "The \"standard deviation\" parameter is to "
                    "high relative to the \"average diameter\" "
                    "parameter. This would results with frequent negative "
                    "diameter values. To fix this problem, please change "
                    "the value of those two parameters or define a "
                    "\"minimum diameter cutoff\" bigger than \"0.\""));
    }
  else
    this->dia_min_cutoff = min_cutoff;

  if (max_cutoff < 0.)
    // approx 99% of all diameters are bigger
    this->dia_max_cutoff = std::exp(mu_ln + 2.5 * sigma_ln);
  else
    this->dia_max_cutoff = max_cutoff;

  AssertThrow(this->dia_min_cutoff < this->dia_max_cutoff,
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
      if (temp_diameter > this->dia_min_cutoff &&
          temp_diameter < this->dia_max_cutoff)
        {
          n_created_diameter++;
          this->particle_sizes.push_back(temp_diameter);
        }
    }
}

double
LogNormalDistribution::find_min_diameter()
{
  return this->dia_min_cutoff;
}

double
LogNormalDistribution::find_max_diameter()
{
  return this->dia_max_cutoff;
}

void
LogNormalDistribution::print_psd_declaration_string(
  const unsigned int        particle_type,
  const ConditionalOStream &pcout)
{
  const double z_min = (std::log(this->dia_min_cutoff) - mu_ln) / sigma_ln;
  const double z_max = (std::log(this->dia_max_cutoff) - mu_ln) / sigma_ln;

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
  : Distribution(
      Parameters::Lagrangian::DistributionWeightingType::number_based)
  , diameter_value(d_values)
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
  const unsigned int        &prn_seed,
  // const double               min_cutoff,
  // const double               max_cutoff,
  const Parameters::Lagrangian::DistributionWeightingType
    &distribution_weighting_type)
  : Distribution(distribution_weighting_type)
  , diameter_custom_values(d_list)
  , gen(prn_seed)
{
  AssertThrow(distribution_weighting_type ==
                DistributionWeightingType::volume_based,
              ExcMessage(
                "Custom distribution only supports volume-based "
                "distribution. Please set the 'distribution weighting basis' "
                "parameter to `volume`."));

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
