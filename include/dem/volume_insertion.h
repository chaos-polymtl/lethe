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
 */

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#ifndef volume_insertion_h
#  define volume_insertion_h

template <int dim>
class VolumeInsertion : public Insertion<dim>
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
  VolumeInsertion(const DEMSolverParameters<dim> &dem_parameters,
                  const double                    maximum_particle_diameter);

  /**
   * Carries out the volume insertion of particles.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   *
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
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
                                 const double         random_number_range,
                                 const int            random_number_seed);

  /**
   * Converts id of particles to volume insertion location
   *
   * @param insertion_location Insertion location of the particle
   * @param id Particle_id
   * @param random_number1 A random number to create randomness in volume insertion
   * @param random_number2 A random number to create randomness in volume insertion
   * @param insertion_information DEM insertion parameters declared in the .prm
   * file
   */
  void
  find_insertion_location_volume(
    Point<dim>                                  &insertion_location,
    const unsigned int                           id,
    const double                                 random_number1,
    const double                                 random_number2,
    const Parameters::Lagrangian::InsertionInfo &insertion_information);

  unsigned int current_inserting_particle_type;

  // Number of particles of each type that remain to be inserted in the
  // upcoming insertion steps
  unsigned int particles_of_each_type_remaining;
};
#endif /* volume_insertion_h */
