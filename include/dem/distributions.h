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
#include <unordered_map>

#ifndef distributions_h
#  define distributions_h

class Distribution
{
public:
  std::vector<double> particle_sizes;

  /**
   * Carries out the size sampling of particles. This is the base class of
   * normal_distribution, log_normal_distribution and list_distribution
   * classes.
   * @param particle_number Number of particle inserted at a given insertion time step.
   * @param particle_type The type of particle being inserted.
   */
  virtual void
  particle_size_sampling(const unsigned int particle_number,
                         const unsigned int particle_type) = 0;
};

class NormalDistribution : public Distribution
{
public:
  /**
   * The constructor stores the parameters necessary to define the normal
   * distribution
   *
   * @param d_averages Average diameters for each type of particle.
   * @param d_standard_deviations Standard deviation of the diameter for each type of particle.
   */
  NormalDistribution(
    const std::unordered_map<unsigned int, double> d_averages,
    const std::unordered_map<unsigned int, double> d_standard_deviations);

  /**
   * Carries out the size sampling of each particle inserted at a insertion
   * time step.
   *
   * @param particle_number Number of particle inserted at a given insertion time step.
   * @param particle_type The type of particle being inserted.
   */
  void
  particle_size_sampling(const unsigned int particle_number,
                         const unsigned int particle_type) override;

private:
  // Average diameters for each particle type.
  std::unordered_map<unsigned int, double> diameter_averages;
  // Standard deviation of the diameter for each particle type.
  std::unordered_map<unsigned int, double> standard_deviations;
};

class UniformDistribution : public Distribution
{
public:
  /**
   * The constructor store the parameters necessary to define the uniform
   * distribution
   *
   * @param d_values Diameter values for each particle type.
   */
  UniformDistribution(const std::unordered_map<unsigned int, double> d_values);

  /**
   * Carries out the size sampling of every particles inserted at a insertion
   * time step.
   *
   * @param particle_number Number of particle inserted at a given insertion time step.
   * @param particle_type The type of particle being inserted.
   */
  void
  particle_size_sampling(const unsigned int particle_number,
                         const unsigned int particle_type) override;

private:
  // Diameter values for each particle type.
  std::unordered_map<unsigned int, double> diameter_values;
};

#endif /* distributions_h */
