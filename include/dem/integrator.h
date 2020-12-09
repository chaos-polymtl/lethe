/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2020 by the Lethe authors
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


#ifndef integration_h
#define integration_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor.h>

#include <deal.II/particles/particle_handler.h>

#include <dem/dem_solver_parameters.h>

using namespace dealii;

/**
 * Base interface for classes that carry out the integration of the velocity and
 * position of particles with inertia
 */

template <int dim>
class Integrator
{
public:
  /**
   * The constructor to the integrator class is currently blank
   * Eventually it might be a good idea to have the integration class contain
   * the index to the velocity property, the force property and the acceleration
   * property manually
   */
  Integrator<dim>()
  {}

  virtual ~Integrator()
  {}

  /**
   * Carries out the integration of the motion of all particles by using
   * the acceleration that is internally stored within the property handler
   * This virtual function that serves as a template for all integrator
   * functions.
   *
   * @param particle_handler The particle handler whose particle motion we wish
   * to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   */
  virtual void
  integrate(Particles::ParticleHandler<dim> &particle_handler,
            Tensor<1, dim>                   body_force,
            double                           time_step) = 0;
};

#endif /* integration_h */
