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
#include <random>

#ifndef distribution_h
#  define distribution_h

class Distribution
{
private:
  /**
   * Carries out the size sampling of particles. This is the base class of
   * normal_distribution, log_normal_distribution and list_distribution
   * classes.
   *
   */
  virtual void
  particle_size_sampling(double particle_number) = 0;

  // Attribute
protected:
  std::vector<double> particle_sizes;
};

#endif /* distribution_h */
