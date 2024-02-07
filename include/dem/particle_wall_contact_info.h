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
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

#ifndef particle_wall_contact_info_h
#  define particle_wall_contact_info_h

/**
 * This struct handles the information related to the calculation of the
 * particle-wall contact force
 */

using namespace dealii;

template <int dim>
class particle_wall_contact_info
{
public:
  /** @brief Construct the particle_wall_contact_info using the particle iterator, the normal vector of the wall,
   *   a point on that boundary and the boundary id. This is the commonly used
   * constructor since it houses all the information required to perform the
   * contact calculation.
   *
   * @param particle The iterator to the particle in contact with the wall
   * @param normal_vector The outward pointing normal vector on the wall
   * @param point_on_boundary A point that lies on the face
   * @param boundary_id The boundary id. This id corresponds to the number attributed to the boundary condition.
   *
   */

  particle_wall_contact_info(const Particles::ParticleIterator<dim> &particle,
                             const Tensor<1, 3>       &normal_vector,
                             const Point<3>           &point_on_boundary,
                             const types::boundary_id &boundary_id)
    : particle(particle)
    , normal_vector(normal_vector)
    , point_on_boundary(point_on_boundary)
    , boundary_id(boundary_id)
    , normal_overlap(0)
    , normal_relative_velocity(0)
    , tangential_overlap({0, 0, 0})
    , tangential_relative_velocity({0, 0, 0})
  {}

  /** @brief Construct the particle_wall_contact_info using only the particle_iterator. This constructor is only used
   *  at one location right now when using floating mesh. It should be
   deprecated in a nearby future
   *
   * @param particle The iterator to the particle in contact with the wall
   *
   */
  particle_wall_contact_info(const Particles::ParticleIterator<dim> &particle)
    : particle(particle)
    , normal_vector({0, 0, 0})
    , point_on_boundary({0, 0, 0})
    , boundary_id(0)
    , normal_overlap(0)
    , normal_relative_velocity(0)
    , tangential_overlap({0, 0, 0})
    , tangential_relative_velocity({0, 0, 0})
  {}



  Particles::ParticleIterator<dim> particle;
  Tensor<1, 3>                     normal_vector;
  Point<3>                         point_on_boundary;
  types::boundary_id               boundary_id;
  double                           normal_overlap;
  double                           normal_relative_velocity;
  Tensor<1, 3>                     tangential_overlap;
  Tensor<1, 3>                     tangential_relative_velocity;
};

#endif /* particle_wall_contact_info_struct_h */
