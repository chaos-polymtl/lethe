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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */
#include <deal.II/base/point.h>

#ifndef INSERTIONINFOSTRUCT_H_
#define INSERTIONINFOSTRUCT_H_

/**
 * This struct handles the information related to the insertion
 * of particles
 */

using namespace dealii;

template <int dim, int spacedim> struct InsertionInfoStruct {
  // Number of insertion phase time steps
  int insertion_steps_number;

  // Number of particles inserted at current step
  int inserted_number_at_step;

  // Frequency of insertion steps
  int insertion_frequency;

  // Ratio of the distance between the centers of two adjacent particles to the
  // particle diameter (for uniform size distribution of particles)
  double distance_threshold;

  // Dimensions of the insertion box
  double x_min;
  double y_min;
  double z_min;
  double x_max;
  double y_max;
  double z_max;
};

#endif /* INSERTIONINFOSTRUCT_H_ */
