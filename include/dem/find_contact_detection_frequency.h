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

#include <deal.II/particles/particle_handler.h>

#include <dem/dem_properties.h>

#include <vector>

using namespace dealii;

#ifndef find_contact_detection_frequency_h
#  define find_contact_detection_frequency_h

/**
 * Carries out
 *
 * @param
 */

template <int dim>
bool
find_contact_detection_frequency(
  Particles::ParticleHandler<dim> &particle_handler,
  const double &                   neighborhood_threshold,
  const double &                   minimum_cell_size,
  const double &                   dt,
  const double &                   maximum_particle_diameter,
  const double &                   dynamic_contact_search_factor);

#endif
