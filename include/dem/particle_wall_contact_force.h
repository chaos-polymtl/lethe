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

#include <math.h>

#include "dem/dem_solver_parameters.h"
#include "dem/physical_info_struct.h"
#include <iostream>
#include <tuple>
#include <vector>

using namespace dealii;

#ifndef PWCONTACTFORCE_H_
#define PWCONTACTFORCE_H_

template <int dim, int spacedim> class ParticleWallContactForce {
public:
  ParticleWallContactForce<dim, spacedim>();
  void pwLinearCF(
      std::vector<std::tuple<
          std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
          Point<dim>, Point<dim>, double, double, double, Point<dim>, double>>,
      physical_info_struct<dim> &);

  void pwNonLinearCF(
      std::vector<std::tuple<
          std::pair<typename Particles::ParticleIterator<dim, spacedim>, int>,
          Point<dim>, Point<dim>, double, double, double, Point<dim>, double>>,
      physical_info_struct<dim> &);

private:
  int sgn(float);
};

#endif /* PWCONTACTFORCE_H_ */
