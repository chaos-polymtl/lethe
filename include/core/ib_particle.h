/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 -  by the Lethe authors
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
 * Author: Lucka Barbeau, Bruno Blais, Polytechnique Montreal, 2019 -
 */

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>


#ifndef lethe_ib_particle_h
#  define lethe_ib_particle_h

using namespace dealii;

template <int dim>
class IBParticle
{
public:
  Point<dim> position;
  // Translational velocity
  Tensor<1, dim> velocity;
  // Angular velocity
  Tensor<1, 3> omega;

  double radius;

  // Pressure imposition location
  Point<dim> pressure_location;
};

#endif
