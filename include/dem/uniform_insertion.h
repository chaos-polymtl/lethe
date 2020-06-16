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

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#ifndef UNIFORMINSERTION_H_
#  define UNIFORMINSERTION_H_

/**
 * Uniform insertion of particles in a rectangular box
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class UniformInsertion : public Insertion<dim>
{
public:
  UniformInsertion<dim>(const DEMSolverParameters<dim> &dem_parameters);

  // UPDATE
  /**
   * Carries out the insertion of particles by discretizing and looping over the
   * insertion box, finding the initial position of particles in the insertion
   * box, their id and other initial properties of the particles. This virtual
   * function that serves as a template for all insertion functions.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param property_pool property pool of particles
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  insert(Particles::ParticleHandler<dim> &                particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;

private:
  // add discription
  virtual std::vector<Point<dim>>
  assign_insertion_points(
    const DEMSolverParameters<dim> &dem_parameters) override;

  // Number of remained particles that should be inserted in the upcoming
  // insertion steps
  unsigned int remained_particles;

  // Number of particles that is going to be inserted at each insetion step.This
  // value can change in the last insertion step to reach the desired number of
  // particles
  unsigned int inserted_this_step;

  //  Number of insertion points in the x, y and z directions, respectively
  unsigned int number_of_particles_x_direction, number_of_particles_y_direction,
    number_of_particles_z_direction;
};

#endif /* UNIFORMINSERTION_H_ */
