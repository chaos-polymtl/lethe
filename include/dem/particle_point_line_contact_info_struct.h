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
 *
 */

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_iterator.h>

#ifndef particle_point_line_contact_info_struct_h
#  define particle_point_line_contact_info_struct_h

/**
 * This struct handles the information related to the calculation of the
 * particle-point and particle-line contact forces
 */

using namespace dealii;

template <int dim>
struct particle_point_line_contact_info_struct
{
  Particles::ParticleIterator<dim> particle;
  Point<3>                         point_one;
  Point<3>                         point_two;
};

#endif /* particle_point_line_contact_info_struct_h */
