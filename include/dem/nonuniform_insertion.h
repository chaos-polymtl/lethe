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

#include "dem/dem_solver_parameters.h"
#include "dem/insertion.h"

#ifndef NONUNIFORMINSERTION_H_
#define NONUNIFORMINSERTION_H_

/**
 * Non-uniform insertion of particles in a rectangular box
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim, int spacedim = dim>
class NonUniformInsertion : public Insertion<dim, spacedim> {
public:
  /**
   * The constructor to the insertion class calculates the maximum number of
   * particles (maximum_particle_number) which can be inserted in the insertion
   * box and prints a warning if the number of particles to be inserted at this
   * step exceeds this calculated number Note that in the current version only
   * the insertion in a rectangular box is defined
   *
   * @param dp Particle diameter
   * @param insertion_info_struct Information related to the insertion of
   * particles from the parameter handler file)
   */
  NonUniformInsertion<dim, spacedim>(
      double dp, InsertionInfoStruct<dim, spacedim> insertion_info_struct);

  /**
   * Carries out the insertion of particles by discretizing and looping over the
   * insertion box, finding the initial position of particles in the insertion
   * box, their id and other initial properties of the particles. This virtual
   * function that serves as a template for all insertion functions.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param tr Triangulation to access the cells in which the particles are
   * inserted
   * @param property_pool property pool of particles
   * @param dp Particle diameter
   * @param rhop Density of particles
   * @param insertion_info_struct Information related to the insertion of
   * particles from the parameter handler file)
   */
  virtual void
  insert(Particles::ParticleHandler<dim, spacedim> &particle_handler,
         const Triangulation<dim, spacedim> &tr,
         Particles::PropertyPool &property_pool, double dp, int rhop,
         InsertionInfoStruct<dim, spacedim> insertion_info_struct) override;
};

#endif /* NONUNIFORMINSERTION_H_ */
