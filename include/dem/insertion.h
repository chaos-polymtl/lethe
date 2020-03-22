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

#include <deal.II/base/array_view.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;

#ifndef INSERTION_H_
#define INSERTION_H_

/**
 * Base interface for classes that carry out the insertion of particles in the
 * system
 */

template <int dim> class Insertion {
public:
  Insertion() {}

  virtual ~Insertion() {}

  /**
   * Carries out the insertion of particles by discretizing and looping over the
   * insertion region, finding the initial position of particles in the
   * insertion region, their id and other initial properties of the particles.
   * This virtual class serves as a template for all insertion functions.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param property_pool Property pool of particles
   * @param dem_parameters DEM parameters declared in the .prm file)
   */
  virtual void insert(Particles::ParticleHandler<dim> &particle_handler,
                      const Triangulation<dim> &triangulation,
                      Particles::PropertyPool &property_pool,
                      const DEMSolverParameters<dim> &dem_parameters) = 0;
};

#endif /* INSERTION_H_ */
