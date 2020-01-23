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
 * Author: Bruno Blais, Polytechnique Montreal, 2019
 */

#include <deal.II/particles/particle_handler.h>

#include <dem/integrator.h>
#include <dem/parameters_dem.h>

using namespace dealii;

#ifndef EXPLICITEULERINTEGRATOR_H_
#  define EXPLICITEULERINTEGRATOR_H_


/**
 * Implementation of a classical explicit euler scheme for the integration
 * of the particle motion
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim, int spacedim = dim>
class ExplicitEulerIntegrator : public Integrator<dim, spacedim>
{
public:
  ExplicitEulerIntegrator()
  {}


  /**
   * Carries out the integration of the motion of all particles by using
   * the acceleration with the explicit Euler method
   *
   * @param particle_handler The particle handler whose particle motion we wish to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param time_step The value of the time step used for the integration
   */
  virtual void
  integrate(Particles::ParticleHandler<dim, spacedim> &particle_handler,
            Tensor<1, dim>                             body_force,
            double                                     time_step) override;
};

#endif
