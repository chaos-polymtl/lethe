/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2023 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 *
 */
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#ifndef distributions_h
#  define distributions_h

class Distribution
{
public:
  std::vector<double> particle_sizes;

  /**
   * @brief Carries out the size sampling of particles. This is the base class of
   * NormalDistribution, UniformDistribution and CustomDistribution classes.
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  virtual void
  particle_size_sampling(const unsigned int &particle_number) = 0;

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
};

class NormalDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the normal
   * distribution.
   *
   * @param d_average Average diameters for a certain normal distribution.
   * @param d_standard_deviation Standard deviation of the diameter for a certain
   * normal distribution.
   */
  NormalDistribution(const double &d_average,
                     const double &d_standard_deviation);

  /**
   * @brief Carries out the size sampling of each particle inserted at an insertion
   * time step for the normal distribution.
   *
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &particle_number) override;

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

private:
  /**
   * Average diameter of the normal distribution.
   */
  const double diameter_average;

  /**
   * Standard deviation of distribution of the normal distribution.
   */
  const double standard_deviation;
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
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  void
  particle_size_sampling(const unsigned int &particle_number) override;

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

private:
  /**
   *  The diameter value of the distribution.
   */
  const double diameter_value;
};

class CustomDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the histogram
   * distribution.
   * @param d_list Vector of diameter values.
   * @param d_probabilities Vector of probability values based on volume fraction
   * with respect to each diameter value.
   */
  CustomDistribution(const std::vector<double> &d_list,
                     const std::vector<double> &d_probabilities);

  /**
   * @brief Carries out the size sampling of each particle inserted at an insertion
   * time step for the histogram distribution.
   *
   * @param particle_number Number of particles inserted at a given insertion time
   * step.
   */
  void
  particle_size_sampling(const unsigned int &particle_number) override;

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

private:
  /**
   * Vector containing all the diameters values.
   */
  const std::vector<double> diameter_custom_values;

  /**
   * Vector containing cumulative probabilities associated with de
   * diameter_custom_values vector. The probabilities are based on the volume
   * fraction, not the number of particles.
   */
  std::vector<double> diameter_custom_cumu_prob;
};

#endif /* distributions_h */
