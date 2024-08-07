/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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
 */

#ifndef lethe_particle_particle_contact_info_struct_h
#define lethe_particle_particle_contact_info_struct_h

#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

using namespace dealii;

/**
 * @brief Handle the information related to the calculation of the
 * particle-particle contact force. Notably it is responsible for storing
 * information that has to be preserved over multiple iterations of a contact,
 * namely everything related to tangential overlaps
 */
template <int dim>
struct particle_particle_contact_info
{
  Particles::ParticleIterator<dim> particle_one;
  Particles::ParticleIterator<dim> particle_two;
  Tensor<1, 3>                     tangential_overlap;
};

#endif
