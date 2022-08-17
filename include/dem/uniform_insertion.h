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

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#ifndef uniform_insertion_h
#  define uniform_insertion_h

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
  UniformInsertion<dim>(const DEMSolverParameters<dim> &dem_parameters,
                        const double maximum_particle_diameter);

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
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  insert(Particles::ParticleHandler<dim> &                particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;

private:
  /**
   * Converts id of particles to uniform insertion location
   *
   * @param insertion_location Insertion location of the particle
   * @param id Particle_id
   * @param insertion_information DEM insertion parameters declared in the .prm
   * file
   */
  void
  find_insertion_location_uniform(
    Point<dim> &                                 insertion_location,
    const unsigned int &                         id,
    const Parameters::Lagrangian::InsertionInfo &insertion_information);

  // Number of remained particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remained_particles_of_each_type;

  unsigned int current_inserting_particle_type;
};

#endif /* uniform_insertion_h */
