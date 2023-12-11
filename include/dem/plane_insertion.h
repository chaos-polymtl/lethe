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


#ifndef plane_insertion_h
#  define plane_insertion_h

/**
 * Insertion of particles using cells cut by a plane
 *
 * @note
 *
 * Particle insertion using cells cut by a plane. Locally owned cells that are
 * cut by the plane are flagged. From those flagged cells, we insert a particle
 * at their centroids if they are individually empty (contain no particle). This
 * way, no significant overlap is occurring on the insertion of new particle
 * which can occur with other insertion method.
 */

template <int dim>
class PlaneInsertion : public Insertion<dim>
{
public:
  /**
   * The plane insertion class insert particles using a plane
   * define by a point an a normal vector. This method of insertion can be
   * useful when dealing with a domain close to be fully filled with particle.
   * In this situation, other insertion method have a high risk to create a big
   * overlap between particles on insertion. The plane insertion method mitigate
   * this risk by insertion at the center of empty cells.
   *
   * @param dem_parameters DEM parameters declared in the .prm file
   * @param triangulation Triangulation object used in the simulation.
   */
  PlaneInsertion(const DEMSolverParameters<dim> &dem_parameters,
                 const parallel::distributed::Triangulation<dim> &triangulation,
                 const std::vector<std::shared_ptr<Distribution>>
                   &distribution_object_container);

  /**
   * Carries out the insertion of particles.
   * @param particle_handler The particle handler of particles which are being
   * inserted.
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted.
   * @param dem_parameters DEM parameters declared in the .prm file.
   */
  virtual void
  insert(Particles::ParticleHandler<dim>                 &particle_handler,
         const parallel::distributed::Triangulation<dim> &triangulation,
         const DEMSolverParameters<dim> &dem_parameters) override;

  virtual void
  serialize(boost::archive::text_oarchive &ar, const unsigned int) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

  virtual void
  deserialize(boost::archive::text_iarchive &ar, const unsigned int) override
  {
    ar &particles_of_each_type_remaining &current_inserting_particle_type;
  }

private:
  /**
   * @brief Store the cells that are cut by the plane.
   *
   * @param triangulation Triangulation to access the cells in which the
   * particles are inserted
   * @param plane_point Point which define the plane
   * @param plane_normal_vector Vector which define the normal direction of
   * the plane.
   */
  void
  find_inplane_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    Point<3>                                         plane_point,
    Tensor<1, 3>                                     plane_normal_vector);

  /**
   * @brief Store the location of the centers of all the cells that are in the plane
   */
  void
  find_centers_of_inplane_cells();


  std::set<typename Triangulation<dim>::active_cell_iterator>
                                               plane_cells_for_insertion;
  double                                       maximum_range_for_randomness;
  int                                          particles_of_each_type_remaining;
  unsigned int                                 current_inserting_particle_type;
  std::unordered_map<unsigned int, Point<dim>> cells_centers;
  bool                                         mark_for_update;
  boost::signals2::connection                  change_to_triangulation;
};


#endif /* plane_insertion_h */
