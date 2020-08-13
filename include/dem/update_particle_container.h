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

#include <deal.II/particles/particle_handler.h>

using namespace dealii;

#ifndef UPDATEPARTICLECONTAINER_H_
#  define UPDATEPARTICLECONTAINER_H_

/**
 * Updates the iterators to local particles in a map of particles
 * (particle_container)
 *
 * @param particle_container A map of particles which is used to update
 * the iterators to particles in pp and pw fine search outputs after calling
 * sort particles into cells function
 * @param particle_handler Particle handler to access all the particles in the
 * system
 */

template <int dim>
void
update_particle_container(
  std::map<int, Particles::ParticleIterator<dim>> &particle_container,
  const Particles::ParticleHandler<dim> *          particle_handler)
{
  particle_container.clear();

  for (auto particle_iterator = particle_handler->begin();
       particle_iterator != particle_handler->end();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }

  for (auto particle_iterator = particle_handler->begin_ghost();
       particle_iterator != particle_handler->end_ghost();
       ++particle_iterator)
    {
      particle_container[particle_iterator->get_id()] = particle_iterator;
    }
}

#endif /* UPDATEPARTICLECONTAINER_H_ */
