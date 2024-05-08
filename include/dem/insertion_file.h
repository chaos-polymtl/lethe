/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/insertion.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#ifndef lethe_insertion_file_h
#  define lethe_insertion_file_h

template <int dim>
class InsertionFile : public Insertion<dim>
{
public:
  InsertionFile(const DEMSolverParameters<dim> &dem_parameters,
                const std::vector<std::shared_ptr<Distribution>>
                  &distribution_object_container);

  /**
   * @brief The InsertionFile class inserts particles using data stored in a file.
   * This allows the insertion of any number of particles at a well-controlled
   * location with any diameter value, translation and angular velocity.
   *
   * @param particle_handler The particle handler of particles which are being
   * inserted
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;


  /**
   * @brief Carries out assigning the properties of inserted particles specifically
   * for the file insertion method. In this method, the initial translation
   * and angular velocities and the diameter of each particles is set.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   *
   * @param inserted_this_step_this_proc Number of particles that are inserted
   * at each insertion step on each processor. This value can change in the last
   * insertion step to reach the desired number of particles
   *
   * @param particles_data Contains the particles
   *
   * @param particle_properties Properties of all inserted particles at each insertion step
   */
  void
  assign_particle_properties_for_file_insertion(
    const DEMSolverParameters<dim>             &dem_parameters,
    const unsigned int                         &inserted_this_step_this_proc,
    std::map<std::string, std::vector<double>> &particles_data,
    std::vector<std::vector<double>>           &particle_properties);



  /**
   * @brief Serialize the list insertion object to an output archive.
   *
   * @param ar Output archive where the attributes are stored.
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar, const unsigned int) override
  {
    ar &remaining_particles_of_each_type &current_inserting_particle_type;
  }

  /**
   * @brief Deserialize an input archive to the list insertion object.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar, const unsigned int) override
  {
    ar &remaining_particles_of_each_type &current_inserting_particle_type;
  }

  // Number of remaining particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remaining_particles_of_each_type;
  unsigned int current_inserting_particle_type;

  // File name where the particles properties are stored.
  std::string file_name;
};

#endif /* lethe_insertion_file_h */
