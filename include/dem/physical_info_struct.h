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

#include <deal.II/base/tensor.h>

#ifndef PHYSICALINFOSTRUCT_H_
#define PHYSICALINFOSTRUCT_H_

/**
 * This struct handles the information related to the physical
 * properties of the simulation including particle and wall properties
 */

using namespace dealii;

template <int dim> struct physical_info_struct {
  double particle_diameter;
  int particle_density;
  int Young_modulus_particle;
  int Young_modulus_wall;
  double Poisson_ratio_particle;
  double Poisson_ratio_wall;
  double restitution_coefficient_particle;
  double restitution_coefficient_wall;
  double friction_coefficient_particle;
  double friction_coefficient_wall;
  double rolling_friction_coefficient_particle;
  double rolling_friction_coefficient_wall;
};

#endif /* PHYSICALINFOSTRUCT_H_ */
