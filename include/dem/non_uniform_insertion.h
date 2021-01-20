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

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#ifndef nonuniform_insertion_h
#  define nonuniform_insertion_h

/**
 * Non-uniform insertion of particles in a rectangular box
 *
 * @note
 *
 * @author Shahab Golshan, Bruno Blais, Polytechnique Montreal 2019-
 */

template <int dim>
class NonUniformInsertion : public Insertion<dim>
{
public:
  /**
   * The constructor investigates if the insertion box is large enough to handle
   * to insert the desired number of particles with the specified insertion
   * parameters. If the insertion box is not adequately large, the number of
   * inserted particles at each insertion step is updated. It also finds the
   * insertion points in each direction (number_of_particles_x_direction,
   * number_of_particles_y_direction and number_of_particles_z_direction).
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param maximum_particle_diameter Maximum particle diameter based on values
   * defined in the parameter handler
   */
  NonUniformInsertion<dim>(const DEMSolverParameters<dim> &dem_parameters,
                           const double &maximum_particle_diameter);

  /**
   * Carries out the non-uniform insertion of particles.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   *
   */
  virtual void
  insert(Particles::ParticleHandler<dim> &                particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;

private:
  /**
   * Creates a vector of random numbers with size of particles which are going
   * to be inserted at each insertion step
   *
   * @param random_container A vector of random numbers
   * @param random_number_range The range in which the random numbers will be
   * generated
   * @param random_number_seed random number seed
   */
  void
  create_random_number_container(std::vector<double> &random_container,
                                 const double &       random_number_range,
                                 const int &          random_number_seed);

  /**
   * Creates a vector of insertion points for non-uniform insertion. The output
   * of this function is used as input argument in insert_global_particles
   *
   * @param insertion_points A vector containing insertion locations of particles
   * @param insertion_information DEM insertion parameters declared in the .prm
   * file
   * @param communicator MPI communicator
   */
  virtual void
  assign_insertion_points(
    std::vector<Point<dim>> &                    insertion_points,
    const Parameters::Lagrangian::InsertionInfo &insertion_information,
    const MPI_Comm &                             communicator) override;

  unsigned int current_inserting_particle_type;

  // Number of particles of each type that remain to be inserted in the
  // upcoming insertion steps
  unsigned int remained_particles_of_each_type;
};

#endif /* nonuniform_insertion_h */
