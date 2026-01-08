// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_distributions_h
#define lethe_distributions_h

#include <core/parameters_lagrangian.h>

#include <random>

class Distribution
{
public:
  std::vector<double> particle_sizes;

  /**
   * @brief Default destructor.
   */
  virtual ~Distribution() = default;

  /**
   * @brief Carries out the size sampling of particles. This is the base class of
   * NormalDistribution, UniformDistribution and CustomDistribution classes.
   * @param number_of_particles Number of particle inserted at a given insertion time step.
   */
  virtual void
  particle_size_sampling(const unsigned int &number_of_particles) = 0;

  /**
   * @brief Return the minimum diameter for a certain distribution.
   *
   * @return The minimum diameter of a certain distribution.
   */
  virtual double
  find_min_diameter() = 0;

  /**
   * @brief Return the maximum diameter for a certain distribution.
   *
   * @return The maximum diameter of a certain distribution.
   */
  virtual double
  find_max_diameter() = 0;

  /**
   * @brief Print the declaration string relative to the particle size
   * distribution used.
   *
   * @param particle_type Particle type of the distribution
   * @param pcout Parallel conditional output stream used to print the
   * information
   */
  virtual void
  print_psd_declaration_string(const unsigned int        particle_type,
                               const ConditionalOStream &pcout) = 0;
};

class NormalDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the normal
   * distribution.
   *
   * @param[in] d_average Average diameters for a certain normal distribution.
   * @param[in] d_standard_deviation Standard deviation of the diameter for a
   * certain normal distribution.
   * @param[in] prn_seed Pseudo-random number seed for the diameter generation.
   * @param[in] min_cutoff Minimum cutoff diameter.
   * @param[in] max_cutoff Maximum cutoff diameter.
   */
  NormalDistribution(const double       &d_average,
                     const double       &d_standard_deviation,
                     const unsigned int &prn_seed,
                     const double       &min_cutoff,
                     const double       &max_cutoff);

  /**
   * @brief Carries out the size sampling of each particle inserted at an insertion
   * time step for the normal distribution.
   *
   * @param[in] number_of_particles Number of particle inserted at a given
   * insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &number_of_particles) override;

  /**
   * @brief Find the minimum diameter a normal distribution.
   *
   * @return The minimum diameter of a normal distribution.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Find the maximum diameter of the normal distribution.
   *
   * @return The maximum diameter of the normal distribution.
   */
  double
  find_max_diameter() override;

  /**
   * @brief Print the declaration string relative to the particle size
   * distribution used.
   *
   * @param particle_type Particle type of the distribution
   * @param pcout Parallel conditional output stream used to print the
   * information
   */
  void
  print_psd_declaration_string(const unsigned int        particle_type,
                               const ConditionalOStream &pcout) override;

private:
  /**
   * @brief Average diameter of the normal distribution.
   */
  const double diameter_average;

  /**
   * @brief Standard deviation of distribution of the normal distribution.
   */
  const double standard_deviation;

  /**
   * @brief Random number generator for the diameter selection.
   */
  std::mt19937 gen;

  /**
   * @brief Minimal cut off diameters values.
   */
  double dia_min_cutoff;

  /**
   * @brief Maximum cut off diameters values.
   */
  double dia_max_cutoff;
};


class LogNormalDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the normal
   * distribution.
   *
   * @param[in] d_average Average diameters for a certain normal distribution.
   * @param[in] d_standard_deviation Standard deviation of the diameter for a
   * certain normal distribution.
   * @param[in] prn_seed Pseudo-random number seed for the diameter generation.
   * @param[in] min_cutoff Minimum cutoff diameter.
   * @param[in] max_cutoff Maximum cutoff diameter.
   */
  LogNormalDistribution(const double       &d_average,
                        const double       &d_standard_deviation,
                        const unsigned int &prn_seed,
                        const double        min_cutoff,
                        const double        max_cutoff);

  /**
   * @brief Carries out the size sampling of each particle inserted at an insertion
   * time step for the normal distribution.
   *
   * @param[in] number_of_particles Number of particle inserted at a given
   * insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &number_of_particles) override;

  /**
   * @brief Find the minimum diameter a normal distribution.
   *
   * @return The minimum diameter of a normal distribution.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Find the maximum diameter of the normal distribution.
   *
   * @return The maximum diameter of the normal distribution.
   */
  double
  find_max_diameter() override;

  /**
   * @brief Print the declaration string relative to the particle size
   * distribution used.
   *
   * @param particle_type Particle type of the distribution
   * @param pcout Parallel conditional output stream used to print the
   * information
   */
  void
  print_psd_declaration_string(const unsigned int        particle_type,
                               const ConditionalOStream &pcout) override;

private:
  /**
   * @brief Standard deviation of distribution of the normal distribution.
   */
  const double sigma_ln;

  /**
   * @brief Average diameter of the normal distribution.
   */
  const double mu_ln;

  /**
   * @brief Random number generator for the diameter selection.
   */
  std::mt19937 gen;

  /**
   * @brief Cut off diameters values.
   */
  double dia_min_cutoff, dia_max_cutoff;
};

class UniformDistribution : public Distribution
{
public:
  /**
   * @brief The constructor store the parameters necessary to define the uniform
   * distribution.
   *
   * @param d_values Diameter values for each particle type.
   */
  UniformDistribution(const double &d_values);

  /**
   * @brief Carries out the size sampling of every particles inserted at an insertion
   * time step for the uniform distribution.
   *
   * @param number_of_particles Number of particle inserted at a given insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &number_of_particles) override;

  /**
   * @brief Find the minimum diameter of the uniform distribution.
   *
   * @return The diameter of the distribution.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Find the maximum diameter of the uniform distribution.
   *
   * @return The diameter of the distribution.
   */
  double
  find_max_diameter() override;


  /**
   * @brief Print the declaration string relative to the particle size
   * distribution used.
   *
   * @param particle_type Particle type of the distribution
   * @param pcout Parallel conditional output stream used to print the
   * information
   */
  void
  print_psd_declaration_string(const unsigned int        particle_type,
                               const ConditionalOStream &pcout) override;

private:
  /**
   * @brief The diameter value of the distribution.
   */
  const double diameter_value;
};

class CustomDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the histogram
   * distribution.
   *
   * @param[in] d_list Vector of diameter values.
   * @param[in] d_probabilities Vector of probability values based on volume
   * fraction with respect to each diameter value.
   * @param[in] prn_seed Pseudo-random number seed for the diameter generation.
   */
  CustomDistribution(const std::vector<double> &d_list,
                     const std::vector<double> &d_probabilities,
                     const unsigned int        &prn_seed);

  /**
   * @brief Carries out the size sampling of each particle inserted at an insertion
   * time step for the histogram distribution.
   *
   * @param[in] number_of_particles Number of particles inserted at a given
   * insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &number_of_particles) override;

  /**
   * @brief Find the minimum diameter of the custom distribution.
   *
   * @return The minimum diameter of the distribution.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Find the maximum diameter of the custom distribution.
   *
   * @return The minimum diameter of the distribution.
   */
  double
  find_max_diameter() override;

  /**
   * @brief Print the declaration string relative to the particle size
   * distribution used.
   *
   * @param particle_type Particle type of the distribution
   * @param pcout Parallel conditional output stream used to print the
   * information
   */
  void
  print_psd_declaration_string(const unsigned int        particle_type,
                               const ConditionalOStream &pcout) override;

private:
  /**
   * Vector containing all the diameters values.
   */
  const std::vector<double> diameter_custom_values;

  /**
   * @brief Vector containing cumulative probabilities associated with de
   * diameter_custom_values vector. The probabilities are based on the volume
   * fraction, not the number of particles.
   */
  std::vector<double> diameter_custom_cumul_prob;

  /**
   * @brief Random number generator for the diameter selection.
   */
  std::mt19937 gen;
};

using namespace Parameters::Lagrangian;
inline void
setup_distributions(const LagrangianPhysicalProperties &lpp,
                    std::vector<std::shared_ptr<Distribution>>
                                &size_distribution_object_container,
                    double      &maximum_particle_diameter,
                    unsigned int mpi_process_id,
                    const ConditionalOStream &pcout)
{
  for (unsigned int particle_type = 0; particle_type < lpp.particle_type_number;
       particle_type++)
    {
      switch (lpp.distribution_type.at(particle_type))
        {
          case SizeDistributionType::uniform:
            size_distribution_object_container[particle_type] =
              std::make_shared<UniformDistribution>(
                lpp.particle_average_diameter.at(particle_type));
            break;
          case SizeDistributionType::normal:
            size_distribution_object_container[particle_type] =
              std::make_shared<NormalDistribution>(
                lpp.particle_average_diameter.at(particle_type),
                lpp.particle_size_std.at(particle_type),
                lpp.seed_for_distributions[particle_type] + mpi_process_id,
                lpp.diameter_min_cutoff.at(particle_type),
                lpp.diameter_max_cutoff.at(particle_type));
            break;
          case SizeDistributionType::lognormal:
            size_distribution_object_container[particle_type] =
              std::make_shared<LogNormalDistribution>(
                lpp.particle_average_diameter.at(particle_type),
                lpp.particle_size_std.at(particle_type),
                lpp.seed_for_distributions[particle_type] + mpi_process_id,
                lpp.diameter_min_cutoff.at(particle_type),
                lpp.diameter_max_cutoff.at(particle_type));
            break;
          case SizeDistributionType::custom:
            size_distribution_object_container[particle_type] =
              std::make_shared<CustomDistribution>(
                lpp.particle_custom_diameter.at(particle_type),
                lpp.particle_custom_probability.at(particle_type),
                lpp.seed_for_distributions[particle_type] + mpi_process_id);
            break;
        }
      size_distribution_object_container[particle_type]
        ->print_psd_declaration_string(particle_type, pcout);

      maximum_particle_diameter = std::max(
        maximum_particle_diameter,
        size_distribution_object_container[particle_type]->find_max_diameter());
    }
}

#endif
