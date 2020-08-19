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

#include <deal.II/grid/filtered_iterator.h>
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
#  define INSERTION_H_

/**
 * Base interface for classes that carry out the insertion of particles in the
 * system
 */

template <int dim>
class Insertion
{
public:
  /**
   * Carries out the insertion of particles. This is the base class of
   * uniform_insertion and non_uniform_insertion classes.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  insert(Particles::ParticleHandler<dim> &                particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &                 dem_parameters) = 0;

protected:
  /**
   * Carries out assigning the properties of inserted particles.
   *
   * @param inserted_this_step Number of particles that are inserted
   * at each insetion step.
   * @param remained_particles Number of remained particles going to be inserted
   * in future insertions
   */
  void
  print_insertion_info(const unsigned int inserted_this_step,
                       const unsigned int remained_particles);

  /**
   * Carries out assigning the properties of inserted particles.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param inserted_this_step Number of particles that are inserted
   * at each insetion step. This value can change in the last insertion step to
   * reach the desired number of particles
   * @param inserted_so_far Number of particles are already inserted
   */
  std::vector<std::vector<double>>
  assign_particle_properties(const DEMSolverParameters<dim> &dem_parameters,
                             const unsigned int &            inserted_this_step,
                             const unsigned int &            inserted_so_far);

  /**
   * Carries out finding the insertion points of inserted particles.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual std::vector<Point<dim>>
  assign_insertion_points(const DEMSolverParameters<dim> &dem_parameters) = 0;
};

#endif /* INSERTION_H_ */
