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

#include "dem/parameters_dem.h"
using namespace dealii;

#ifndef INTEGRATION_H_
#  define INTEGRATION_H_

template <int dim, int spacedim = dim>
class Integrator
{
public:
  Integrator();

  virtual ~Integrator(){};

  /**
   * Carries out the integration of the motion of all particles by using
   * the acceleration that is internally stored within the property handler
   * This virtual function that serves as a template for all integrator
   * function.
   *
   * @param particle_handler The particle handler whose particle motion we wish to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   */
  virtual void
  integrate(Particles::ParticleHandler<dim, spacedim> &particle_handler,
            Tensor<1, dim>                             body_force,
            double                                     time_step)
  {}

  // void
  // rk2Integration(Particles::ParticleHandler<dim, spacedim> &particle_handler,
  //               Point<dim>                                 body_force,
  //               float                                      time_step);

  void
  gearIntegration();
};

#endif /* INTEGRATION_H_ */
