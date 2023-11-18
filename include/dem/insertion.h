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
 */

#include <core/dem_properties.h>

#include <dem/dem_solver_parameters.h>
#include <dem/normal_distribution.h>

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

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace dealii;

#ifndef insertion_h
#  define insertion_h

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
   * volume_insertion, plane_insertion and list_insertion classes.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  Insertion(const DEMSolverParameters<dim> &dem_parameters);


  // Destructor
  /**
   * This function is override by volume_insertion, plane_insertion and
   * list_insertion class to insert particles.
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
         const DEMSolverParameters<dim>                  &dem_parameters) = 0;

protected:
  /**
   * Carries out assigning the properties of inserted particles.
   *
   * @param inserted_this_step Number of particles that are inserted
   * at each insertion step.
   * @param remained_particles Number of remained particles going to be inserted
   * in future insertions
   * @param particle_type Inserted particle type
   * @param pcout Printing in parallel
   */
  void
  print_insertion_info(const unsigned int       &inserted_this_step,
                       const unsigned int       &remained_particles,
                       const unsigned int       &particle_type,
                       const ConditionalOStream &pcout);

  /**
   * Carries out assigning the properties of inserted particles.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param inserted_this_step_this_proc Number of particles that are inserted
   * at each insertion step on each processor. This value can change in the last
   * insertion step to reach the desired number of particles
   * @param current_inserting_particle_type Type of inserting particles
   * @param particle_properties Properties of all inserted particles at each insertion step
   */
  void
  assign_particle_properties(
    const DEMSolverParameters<dim>   &dem_parameters,
    const unsigned int               &inserted_this_step_this_proc,
    const unsigned int               &current_inserting_particle_type,
    std::vector<std::vector<double>> &particle_properties);

  /**
   * @brief Carries out finding the maximum number of inserted particles based on the
   * insertion box size. If the requested number of particles for insertion in
   * each insertion step is larger than this maximum, it is limited to this
   * value and a warning is printed.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param pcout Printing in parallel
   */
  void
  calculate_insertion_domain_maximum_particle_number(
    const DEMSolverParameters<dim> &dem_parameters,
    const ConditionalOStream       &pcout);

  // Number of particles that is going to be inserted at each insertion
  // step.This value can change in the last insertion step to reach the desired
  // number of particles
  unsigned int inserted_this_step;

  // Number of insertion points in the x, y and z directions
  std::vector<int> number_of_particles_directions;

  // Minimum and maximum number of inserted particles based on the insertion box
  // size and the direction order (it means that axis 0 is not necessarily x
  // etc...) It depends of the order of the insertion direction.
  std::vector<double> axis_min, axis_max;

  // Maximum particle diameter
  double maximum_diameter;

  // Inserted number of particles at this step on this processor
  unsigned int inserted_this_step_this_proc;

  // A vector of vectors, which contains all the properties of all inserted
  // particles at each insertion step
  std::vector<std::vector<double>> particle_properties;

  // A distribution object that carries out the attribution of diameter to every
  // particle during an insertion time step
  std::shared_ptr<Distribution> distribution_object;

private:
  // Stores particle diameter
  std::vector<double> particle_sizes;
};

#endif /* insertion_h */
