// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_insertion_h
#define lethe_insertion_h

#include <dem/dem_solver_parameters.h>
#include <dem/distributions.h>

#include <deal.II/base/data_out_base.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace dealii;

/**
 * @brief Base interface for classes that carry out the insertion of particles
 * in the system.
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class Insertion
{
public:
  /**
   * @brief Carries out the insertion of particles. This is the base class of
   * volume_insertion, plane_insertion, list_insertion and file_insertion
   * classes.
   *
   * @param size_distribution_object_container Contains all distribution for each
   * particle type
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param dem_parameters DEM parameters declared in the .prm file
   */
  Insertion(const std::vector<std::shared_ptr<Distribution>>
              &size_distribution_object_container,
            const parallel::distributed::Triangulation<dim> &triangulation,
            const DEMSolverParameters<dim>                  &dem_parameters);
  /**
   * @brief Carries out the destruction of the insertion class.
   *
   */
  virtual ~Insertion()
  {
    // This boost signal needs to be disconnected before being destroyed,
    // otherwise a segfault could occur when the triangulation will be
    // destroyed.
    this->change_to_triangulation.disconnect();
  };

  /**
   * @brief This function is overridden by the volume_insertion, plane_insertion,
   * list_insertion and file_insertion classes.
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

  /**
   * @brief Serialize the insertion object to an output archive. Is being used
   * when checkpointing a simulation. This function is overridden by
   * volume_insertion, plane_insertion, list_insertion and file_insertion
   * classes.
   *
   * @param ar Output archive where the attributes are stored.
   *
   */
  virtual void
  serialize(boost::archive::text_oarchive &ar) = 0;


  /**
   * @brief Deserialize an input archive to the insertion object. Is being used
   * when restarting a simulation. This function is overridden by
   * volume_insertion, plane_insertion, list_insertion and from_file_insertion
   * classes.
   *
   * @param ar Input archive where the attributes are stored.
   *
   */
  virtual void
  deserialize(boost::archive::text_iarchive &ar) = 0;

protected:
  /**
   * @brief Print information about the particles that have been inserted during an
   * insertion time step.
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
   * @brief Carries out assigning the properties of inserted particles.
   *
   * @param dem_parameters DEM parameters declared in the prm file.
   * @param inserted_this_step_this_proc Number of particles that are inserted.
   * at each insertion step on each processor. This value can change in the last
   * insertion step to reach the desired number of particles.
   * @param current_inserting_particle_type Type of inserting particles.
   * @param insertion_points Insertion points of all inserted particles at this
   * insertion step.
   * @param particle_properties Properties of all inserted particles at this
   * insertion step.
   */
  void
  assign_particle_properties(
    const DEMSolverParameters<dim>   &dem_parameters,
    const unsigned int               &inserted_this_step_this_proc,
    const unsigned int               &current_inserting_particle_type,
    const std::vector<Point<dim>>    &insertion_points,
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

  /**
   * @brief Find every cell that are completely and partially inside de the
   * clearing box. For the cell completely inside the box, all the particles
   * will be deleted. For those partially inside, a fine verification is
   * required.
   *
   * @param triangulation Triangulation to access the cells in which the
   * particles deleted.
   */
  void
  find_cells_in_removing_box(
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * @brief Remove every particle located inside a predefined zone. Currently,
   * the zone is defined by a box.
   *
   * @param particle_handler The particle handler of particles which are being
   * removed.
   */
  void
  remove_particles_in_box(Particles::ParticleHandler<dim> &particle_handler);

  // Number of particles inserted at each insertion time step. This value can
  // change in the last insertion step to reach the desired number of particles
  unsigned int inserted_this_step;

  // Number of insertion points in the x, y and z directions
  std::vector<int> number_of_particles_directions;

  // Minimum and maximum number of inserted particles based on the insertion box
  // size and the direction order (it means that axis 0 is not necessarily x
  // etc...) It depends on the order of the insertion direction.
  std::vector<double> axis_min, axis_max;

  // Maximum particle diameter
  double maximum_diameter;

  // Inserted number of particles at this step on this processor
  unsigned int inserted_this_step_this_proc;

  // A distribution object that carries out the attribution of diameter to every
  // particle during an insertion time step
  std::vector<std::shared_ptr<Distribution>> distributions_objects;

  // Clearing bool
  const bool removing_particles_in_region;

  // Clearing box coordinates
  Point<dim> p_min, p_max;

  // Cell iterator containers for clearing
  std::set<typename Triangulation<dim>::active_cell_iterator> in_removal_box;
  std::set<typename Triangulation<dim>::active_cell_iterator>
    edge_of_removal_box;

  // For when the triangulation has changed (i.e. when load balancing)
  bool                        mark_for_update;
  boost::signals2::connection change_to_triangulation;

private:
  // Stores particles diameters
  std::vector<double> particle_sizes;
};

#endif /* insertion_h */
