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

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <tuple>
#include <vector>

#include "dem/parameters_dem.h"

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#  define PWCONTACTFORCE_H_

class ParticleWallContactForce
{
public:
  ParticleWallContactForce();
  void pwLinearCF(
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>, int>,
                           Point<3>,
                           Point<3>,
                           double,
                           double,
                           double,
                           Point<3>,
                           double>>,
    Particles::ParticleHandler<3, 3> &,
    ParametersDEM<3>);

private:
  double vecValue(Point<3>);
  int
  sgn(float);
};

#endif /* PWCONTACTFORCE_H_ */
