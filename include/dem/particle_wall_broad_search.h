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

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>
#include <dem/disable_contacts.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>


using namespace dealii;

#ifndef particle_wall_broad_search_h
#  define particle_wall_broad_search_h

/**
 * This class is used for broad particle-wall contact search. Broad search
 * is used to obtain all the particles located at boundary cells
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleWallBroadSearch
{
public:
  ParticleWallBroadSearch<dim>();

  /**
   * Finds unordered map of tuples (tuple of particle located in
   * boundary cells, normal vector of the boundary face, a
   * point on the face and the corresponding boundary cell) which shows the
   * candidate particle-wall collision pairs. These collision candidates will be
   * investigated in the fine search to check if they are in contact or not
   *
   * @param boundary_cells_information Information of the boundary cells and
   * faces. This is the output of the FindBoundaryCellsInformation class
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param particle_wall_contact_candidates A two-layered unordered map of tuples. Each
   * tuple contains a particle located near boundaries, the normal vector of
   * the corresponding face boundary, a point on the boundary and the boundary
   * cell. The contact pair is used in the fine search
   */

  void
  find_particle_wall_contact_pairs(
    const std::map<int, boundary_cells_info_struct<dim>>
      &                                    boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates
      &particle_wall_contact_candidates);

  void
  find_particle_wall_contact_pairs(
    const std::map<int, boundary_cells_info_struct<dim>>
      &                                    boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates
      &                         particle_wall_contact_candidates,
    const DisableContacts<dim> &disable_contacts_object);

  /**
   * Finds a two-layered unordered map of particle iterators which shows the
   * candidate particle-floating wall collision candidates. These collision
   * pairs will be investigated in the fine search to check if they are in
   * contact or not
   *
   * @param boundary_cells_for_floating_walls Boundary cells located adjacent to
   * floating walls
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param floating_wall_properties Properties of the floating walls specified
   * in the parameter handler file
   * @param simulation_time Simulation time
   * @param particle_floating_wall_candidates Output of particle-floating wall broad search
   * which contains all the particle-floating wall collision candidates
   */

  void
  find_particle_floating_wall_contact_pairs(
    const std::unordered_map<
      unsigned int,
      std::set<typename Triangulation<dim>::active_cell_iterator>>
      &                                    boundary_cells_for_floating_walls,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double                                      simulation_time,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
      &particle_floating_wall_candidates);

  void
  find_particle_floating_wall_contact_pairs(
    const std::unordered_map<
      unsigned int,
      std::set<typename Triangulation<dim>::active_cell_iterator>>
      &                                    boundary_cells_for_floating_walls,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double                                      simulation_time,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
      &                         particle_floating_wall_candidates,
    const DisableContacts<dim> &disable_contacts_object);

  /**
   * Finds a two-layered unordered map
   * (particle_floating_mesh_contact_candidates) of particle iterators that
   * shows the candidate particle-floating mesh collision candidates. These
   * collision pairs will be investigated in the fine search to check if they
   * are in contact or not
   *
   * @param floating_mesh_information Information of the floating mesh mapped in the
   * background triangulation
   * @param particle_handler
   * @param particle_floating_mesh_contact_candidates Particle-floating mesh contact
   * candidates
   * @param cells_total_neighbor_list A container in which all the neighbor cells
   * of the local cells are stored
   */

  void
  particle_floating_mesh_contact_search(
    const typename DEM::dem_data_structures<dim>::floating_mesh_information
      &                                    floating_mesh_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
      &particle_floating_mesh_contact_candidates,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
      &cells_total_neighbor_list);

  void
  particle_floating_mesh_contact_search(
    const typename DEM::dem_data_structures<dim>::floating_mesh_information
      &                                    floating_mesh_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
      &particle_floating_mesh_contact_candidates,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
      &                         cells_total_neighbor_list,
    const DisableContacts<dim> &disable_contacts_object);
};

#endif /* particle_wall_broad_search_h */
