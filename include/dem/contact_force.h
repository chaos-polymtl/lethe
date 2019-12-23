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

#ifndef CONTACTFORCE_H_
#  define CONTACTFORCE_H_

class ContactForce
{
public:
  ContactForce();
  void linearCF(
    std::vector<std::tuple<std::pair<Particles::ParticleIterator<3, 3>,
                                     Particles::ParticleIterator<3, 3>>,
                           double,
                           Point<3>,
                           double,
                           Point<3>,
                           double,
                           double>>,
    Particles::ParticleHandler<3, 3> &,
    ParametersDEM<3>);


private:
  int
  sgn(float);
};

#endif /* CONTACTFORCE_H_ */
