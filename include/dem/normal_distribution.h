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

#include <dem/distribution.h>

#ifndef normal_distribution_h
#  define normal_distribution_h

class NormalDistribution : public Distribution
{
  /**
   * The constructor store the parameters necessary to define the normal
   * distribution
   *
   * @param average Average of the normal distribution.
   * @param standard_deviation Stabdard_deviation of the normal distribution.
   * @param particle_number Number of particles to insert at the current
   * insertion time step.
   */
  NormalDistribution(const double normal_distribution_average,
                     const double normal_distribution_standard_deviation);

  /**
   * Carries out the size sampling of every particles inserted at a insertion
   * time step.
   *
   */
  void
  particle_size_sampling(double particle_number) override;

  // Attributes
  double average;
  double standard_deviation;
};
#endif /* normal_distribution_h */
