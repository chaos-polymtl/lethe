/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2021 by the Lethe authors
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

#ifndef lethe_radioactive_particle_h
#define lethe_radioactive_particle_h

/**
 * Contains all properties of a particle at one position.
 */

#include <deal.II/base/point.h>

using namespace dealii;

template <int dim>
class RadioParticle
{
public:
  /**
   * @brief Constructor for the RadioParticle.
   *
   * @param location Position of the particle
   *
   * @param n ID number to the particle
   *
   */
  RadioParticle(Point<dim> &location, int n)
    : position(location)
    , id(n)
  {}

  Point<dim>
  get_position()
  {
    return position;
  }

  Point<dim>
  get_id()
  {
    return id;
  }

private:
  Point<dim> position;
  int        id;
};

#endif // lethe_radioactive_particle_h
