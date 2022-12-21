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

#ifndef particle_wall_contact_info_struct_h
#  define particle_wall_contact_info_struct_h

/**
 * This struct handles the information related to the calculation of the
 * particle-wall contact force
 */

using namespace dealii;

template <int dim>
class particle_wall_contact_info_struct
{
public:
  Particles::ParticleIterator<dim> particle;
  Tensor<1, 3>                     normal_vector;
  Point<3>                         point_on_boundary;
  double                           normal_overlap;
  double                           normal_relative_velocity;
  Tensor<1, 3>                     tangential_overlap;
  Tensor<1, 3>                     tangential_relative_velocity;
  types::boundary_id               boundary_id;
  int                              global_face_id;

  // particle_wall_contact_info_struct(Particles::ParticleIterator<dim>
  // &particle,
  //                                   Tensor<1,3> &normal_vector,
  //                                   Point<3> &point_on_boundary)
};

#endif /* particle_wall_contact_info_struct_h */
