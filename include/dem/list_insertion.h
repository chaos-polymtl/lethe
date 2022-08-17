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

#ifndef lethe_list_insertion_h
#  define lethe_list_insertion_h

/**
 * Insertion of particles using a list of positions
 *
 * @author Bruno Blais, Polytechnique Montreal 2021-
 */

template <int dim>
class ListInsertion : public Insertion<dim>
{
public:
  ListInsertion<dim>(const DEMSolverParameters<dim> &dem_parameters);

  /**
   * @brief The ListInsertion class inserts particles using a list specific position.
   * This allows the insertion of any number of particles at a well-controled
   * location which is especially useful from a testing perspective. The code
   * ensures that the number of positions provided in the x,y (and possibly z)
   * direction is coherent. If more particles than the number of position in the
   * list are requested, the class will continue inserting particles at the
   * insertion frequency using the list of position. There is no mechanism in
   * place that prevents the overlap of these new particles with previous ones.
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

  // Number of remaining particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remaining_particles_of_each_type;
  unsigned int current_inserting_particle_type;

  // Vector of the location used to insert the particles
  std::vector<Point<dim>> insertion_points;
};

#endif /* uniform_insertion_h */
