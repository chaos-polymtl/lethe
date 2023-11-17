/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#include <dem/distribution_base.h>

#include <unordered_map>

#ifndef normal_distribution_h
#  define normal_distribution_h

class NormalDistribution : public Distribution
{
public:
  /**
   * The constructor store the parameters necessary to define the normal
   * distribution
   *
   * @param average Average of the normal distribution.
   * @param standard_deviation Standard_deviation of the normal distribution.
   * @param particle_number Number of particles to insert at the current
   * insertion time step.
   */
  NormalDistribution(
    const std::unordered_map<unsigned int, double> normal_distribution_average,
    const std::unordered_map<unsigned int, double>
      normal_distribution_standard_deviation);

  /**
   * Carries out the size sampling of every particles inserted at a insertion
   * time step.
   *
   */
  void
  particle_size_sampling(const unsigned int particle_number,
                         const unsigned int particle_type) override;

private:
  // Attributes
  std::unordered_map<unsigned int, double> diameter_averages;
  std::unordered_map<unsigned int, double> standard_deviations;
};
#endif /* normal_distribution_h */
