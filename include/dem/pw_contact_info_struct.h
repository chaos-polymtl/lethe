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


#ifndef particle_wall_contact_info_struct_h
#define particle_wall_contact_info_struct_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

/**
 * This struct handles the information related to the calculation of the
 * particle-wall contact force
 */

using namespace dealii;

template <int dim>
struct pw_contact_info_struct
{
  Particles::ParticleIterator<dim> particle;
  Tensor<1, dim>                   normal_vector;
  Point<dim>                       point_on_boundary;
  double                           normal_overlap;
  double                           normal_relative_velocity;
  Tensor<1, dim>                   tangential_overlap;
  Tensor<1, dim>                   tangential_relative_velocity;
  unsigned int                     face_id;
  unsigned int                     boundary_id;
};

#endif /* particle_wall_contact_info_struct_h */
