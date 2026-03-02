// SPDX-FileCopyrightText: Copyright (c) 2021-2022 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

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
