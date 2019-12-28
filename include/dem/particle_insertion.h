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

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#include "dem/parameters_dem.h"


#ifndef PARTICLEINSERTION_H_
#  define PARTICLEINSERTION_H_

template <int dim, int spacedim>
class ParticleInsertion
{
public:
  ParticleInsertion<dim, spacedim>(ParametersDEM<dim>);
  void
  uniformInsertion(dealii::Particles::ParticleHandler<dim, spacedim> &,
                   const dealii::Triangulation<dim, spacedim> &,
                   ParametersDEM<dim>,
                   int &,
                   dealii::Particles::PropertyPool &);
};

#endif /* PARTICLEINSERTION_H_ */
