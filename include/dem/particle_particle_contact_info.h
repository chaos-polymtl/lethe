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
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

#ifndef particle_particle_contact_info_struct_h
#  define particle_particle_contact_info_struct_h

/** @brief
 * This class handles the information related to the calculation of the
 * particle-particle contact force. Notably it is responsible for storing
 * information that has to be preserved over multiple iterations of a contact,
 * namely everything related to tangential overlaps
 */

using namespace dealii;

template <int dim>
class particle_particle_contact_info
{
public:
  /** @brief Constructor for cases where the two particle iterators are known.
   *  This constructor is used everywhere in the regular DEM and unresolved
   * CFD-DEM code
   *
   *  @param particle_one A reference to the iterator of particle one
   *  @param particle_two A reference to the iterator of particle two
   */

  inline particle_particle_contact_info(
    Particles::ParticleIterator<dim> &particle_one,
    Particles::ParticleIterator<dim> &particle_two)
    : particle_one(particle_one)
    , particle_two(particle_two)
  {}

  /** @brief Dummy constructor for cases where the two particle iterators are unknown.
   *  This constructor is only used in the resolved CFD-DEM and should be
   * deprecated eventually.
   */

  inline particle_particle_contact_info()
  {}


  Tensor<1, 3>                     tangential_relative_velocity;
  Tensor<1, 3>                     tangential_overlap;
  Particles::ParticleIterator<dim> particle_one;
  Particles::ParticleIterator<dim> particle_two;
};

#endif /* particle_particle_contact_info_struct_h */
