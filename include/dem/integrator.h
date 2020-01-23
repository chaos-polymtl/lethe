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
  Integrator<dim, spacedim>();

  /**
   * Copy a subset of the rows and columns of another matrix into the current
   * object.
   *
   * @param particle_handler The particle handler whose particle motion we wish to integrate
   * @param body_force A constant volumetric body force applied to all particles
   * @param column_index_set The set of columns of @p matrix from which to
   * extract. @pre The number of elements in @p row_index_set and @p
   * column_index_set shall be equal to the number of rows and columns in the
   * current object. In other words, the current object is not resized for
   * this operation.
   */
  void
  eulerIntegration(Particles::ParticleHandler<dim, spacedim> &particle_handler,
                   Point<dim>                                 body_force,
                   float                                      time_step);
  void
  rk2Integration(Particles::ParticleHandler<dim, spacedim> &particle_handler,
                 Point<dim>                                 body_force,
                 float                                      time_step);
  void
  velocity_verlet_integration(
    Particles::ParticleHandler<dim, spacedim> &particle_handler,
    Point<dim>                                 body_force,
    float                                      time_step);
  void
  gearIntegration();
};

#endif /* INTEGRATION_H_ */
