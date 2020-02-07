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
#include <deal.II/particles/particle_iterator.h>

#ifndef CONTACTINFOSTRUCT_H_
#define CONTACTINFOSTRUCT_H_

using namespace dealii;

template <int dim, int spacedim> struct contact_info_struct {
  double normal_overlap;
  Point<dim> normal_vector;
  double normal_relative_velocity;
  Point<dim> tangential_vector;
  double tangential_relative_velocity;
  double tangential_overlap;
  Particles::ParticleIterator<dim, spacedim> particle_one;
  Particles::ParticleIterator<dim, spacedim> particle_two;
};

#endif /* CONTACTINFOSTRUCT_H_ */
