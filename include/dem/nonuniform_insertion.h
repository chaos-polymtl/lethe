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
   * @param x_min, y_min, z_min, x_max, y_max and z_max Dimensions of the
   * rectangular insertion box
   * @param dp Particle diameter
   * @param inserted_number_at_step Inserted number of particles at each
   * insertion step
   * @param distance_threshold Considered distance between two particles (ratio
   * of distance to diameter of particles)
   */
  NonUniformInsertion<dim, spacedim>(double x_min, double y_min, double z_min,
                                  double x_max, double y_max, double z_max,
                                  double dp, int inserted_number_at_step,
                                  double distance_threshold);

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
   * @param active_particle_number Number of active particles inserted so far
   * @param property_pool property pool of particles
   * @param x_min, y_min, z_min, x_max, y_max and z_max Dimensions of the
   * rectangular insertion box
   * @param dp Particle diameter
   * @param inserted_number_at_step Inserted number of particles at each
   * insertion step
   * @param g gravitational acceleration
   * @param distance_threshold Considered distance between two particles (ratio
   * of distance to diameter of particles)
   */
  virtual void
  insert(Particles::ParticleHandler<dim, spacedim> &particle_handler,
         const Triangulation<dim, spacedim> &tr, int &active_particle_number,
         Particles::PropertyPool &property_pool, double x_min, double y_min,
         double z_min, double x_max, double y_max, double z_max, double dp,
         int inserted_number_at_step, int rhop, Tensor<1, dim> g,
         double distance_threshold) override;
};

#endif /* NONUNIFORMINSERTION_H_ */
