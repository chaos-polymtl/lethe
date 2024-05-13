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
#include <dem/insertion_file.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#ifndef lethe_insertion_clear_and_add_h
#  define lethe_insertion_list_of_file_h

template <int dim>
class InsertionClearAndAdd : public Insertion<dim>
{
public:
  InsertionClearAndAdd(
    const DEMSolverParameters<dim>                  &dem_parameters,
    const parallel::distributed::Triangulation<dim> &triangulation,
    const std::vector<std::shared_ptr<Distribution>>
      &distribution_object_container);

  /**
   * @brief The InsertionClearAndAdd class deletes every particle present in a
   * certain region then inserts particles using data stored in files.
   * This allows the insertion of any number of particles at a
   * well-controlled location with any diameter value, translation and angular
   * velocity. If well used, the insertion is guarantee to not caused an
   * excessive overlap with the particles already present in the triangulation.
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
   * @brief Find every cell that are completely and partially inside de the
   * clearing box. For the cell completely inside the box, all the particles
   * will be deleted. For those partially inside, fine verification is required.
   *
   * @param triangulation Triangulation to access the cells in which the
   * particles deleted.
   */
  virtual void
  find_in_clearing_box_cells(
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * @brief Carries out assigning the properties of inserted particles specifically
   * for the ClearAndAdd insertion method. In this method, the initial
   * translation and angular velocities and the diameter of each particles is
   * set.
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
    ar &remaining_particles_of_each_type &current_inserting_particle_type
      &current_file_id;
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
    ar &remaining_particles_of_each_type &current_inserting_particle_type
      &current_file_id;
  }

  // Number of remaining particles of each type that should be inserted in the
  // upcoming insertion steps
  unsigned int remaining_particles_of_each_type;
  unsigned int current_inserting_particle_type;

  // File id iterator
  unsigned int       current_file_id;
  const unsigned int number_of_files;

  // Files where the particles properties are stored.
  const std::vector<std::string> insertion_files;

  // Clearing box coordinates
  Point<dim> p_min, p_max;

  // Cell iterator containers
  std::set<typename Triangulation<dim>::active_cell_iterator>
    in_the_clearing_box;
  std::set<typename Triangulation<dim>::active_cell_iterator> edge_of_box;

  // For when the triangulation has change (i.e. when load balancing)
  bool                        mark_for_update;
  boost::signals2::connection change_to_triangulation;
};

#endif /* lethe_insertion_clear_and_add_h */
