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

#include <deal.II/particles/particle_handler.h>

#include "read_input_script.h"
using namespace dealii;

#ifndef INTEGRATION_H_
#  define INTEGRATION_H_

class Integration
{
public:
  Integration();
  void eulerIntegration(Particles::ParticleHandler<3, 3> &, ReadInputScript);
  void velVerIntegration(Particles::ParticleHandler<3, 3> &, float);
  void
  gearIntegration();
};

#endif /* INTEGRATION_H_ */
