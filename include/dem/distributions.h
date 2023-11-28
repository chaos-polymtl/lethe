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
   * NormalDistribution, UniformDistribution, LogNormalDistribution and
   * HistogramDistribution classes.
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  virtual void
  particle_size_sampling(const unsigned int particle_number) = 0;

  /**
   * @brief Return the minimum diameter for a certain distribution.
   */
  virtual double
  find_min_diameter() = 0;

  /**
   * @brief Return the maximum diameter for a certain distribution.
   */
  virtual double
  find_max_diameter() = 0;
};

class NormalDistribution : public Distribution
{
public:
  /**
   * @brief The constructor stores the parameters necessary to define the normal
   * distribution
   *
   * @param d_averages Average diameters for each type of particle.
   * @param d_standard_deviations Standard deviation of the diameter for each type of particle.
   */
  NormalDistribution(const double d_averages,
                     const double d_standard_deviations);

  /**
   * @brief Carries out the size sampling of each particle inserted at a insertion
   * time step.
   *
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  void
  particle_size_sampling(const unsigned int particle_number) override;

  /**
   * @brief Return @diameter_values.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Return @diameter_values.
   */
  double
  find_max_diameter() override;

private:
  // Average diameters.
  double diameter_averages;
  // Standard deviation of the diameter.
  double standard_deviations;
};

class UniformDistribution : public Distribution
{
public:
  /**
   * @brief The constructor store the parameters necessary to define the uniform
   * distribution
   *
   * @param d_values Diameter values for each particle type.
   */
  UniformDistribution(const double d_values);

  /**
   * @brief Carries out the size sampling of every particles inserted at a insertion
   * time step.
   *
   * @param particle_number Number of particle inserted at a given insertion time step.
   */
  void
  particle_size_sampling(const unsigned int particle_number) override;


  /**
   * @brief Return @diameter_values.
   */
  double
  find_min_diameter() override;

  /**
   * @brief Return @diameter_values.
   */
  double
  find_max_diameter() override;

private:
  // Diameter values.
  double diameter_values;
};

#endif /* distributions_h */
