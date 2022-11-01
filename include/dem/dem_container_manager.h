/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
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

#include <core/serial_solid.h>

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/localize_contacts.h>
#include <dem/locate_local_particles.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_fine_search.h>

#include <deal.II/particles/particle_handler.h>

using namespace DEM;

#ifndef lethe_dem_container_manager_h
#  define lethe_dem_container_manager_h

template <int dim>
class DEMContainerManager
{
public:
  /**
   * Particle-particle broad search which finds the particle-particle contact
   * pairs candidates, local or ghost, (contact_pair_candidates) which shows the
   * collision pairs. These collision pairs will be used in the fine search
   * to investigate if they are in contact or not.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   */
  void
  execute_particle_particle_broad_search(
    dealii::Particles::ParticleHandler<dim> &particle_handler);

  /**
   * Iterates over a vector of maps (pairs_in_contact) to see if the particles
   * which were in contact in the last time step, are still in contact or not.
   *
   * @param neighborhood_threshold A value which defines the neighbor particles
   */
  void
  execute_particle_particle_fine_search(const double neighborhood_threshold);

  /**
   * @brief Carries out the broad contact detection search using the
   * background triangulation for particle-walls contact
   *
   */
  void
  execute_particle_wall_broad_search(
    const Particles::ParticleHandler<dim> &           particle_handler,
    BoundaryCellsInformation<dim> &                   boundary_cell_object,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const bool has_floating_mesh = false);

  /**
   * @brief Carries out the fine particle-wall contact detection
   *
   */
  void
  execute_particle_wall_fine_search(
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const double                                      neighborhood_threshold,
    const bool has_floating_mesh = false);

  /**
   * Manages to call update_fine_search_candidates() with contact data
   * containers to remove contact repetitions and to add new contact pairs to
   * the contact containers when particles are exchanged between processors.
   * Repeats the update for local particle-particle contacts, ghost
   * particle-particle, particle-wall contacts, particle-floating wall contacts
   * and particle-floating mesh contacts.
   *
   */
  void
  localize_contacts();

  /**
   * Updates the iterators to particles in local-local contact containers. This
   * is essential since sort_particles_into_subdomains_and_cells() and
   * exchange_ghost_particles() functions change the iterator to particles
   * everytime they are called.
   *
   * @tparam dim Dimensionality of the geometry which contains the particles
   * @param particle_handler
   *
   */
  void
  locate_local_particles_in_cells(
    const Particles::ParticleHandler<dim> &particle_handler);

  void
  execute_cell_neighbors_search(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const bool                                       has_floating_mesh = false);

  void
  update_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const bool                                       has_floating_mesh = false)
  {
    cells_local_neighbor_list.clear();
    cells_ghost_neighbor_list.clear();

    execute_cell_neighbors_search(triangulation, has_floating_mesh);
  }

  void
  store_floating_mesh_info(
    const parallel::distributed::Triangulation<dim> &       triangulation,
    std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids);

  inline void
  clear_cells_neighbor_list()
  {
    cells_local_neighbor_list.clear();
    cells_ghost_neighbor_list.clear();
  }

  inline void
  clear_contact_pair_candidates()
  {
    local_contact_pair_candidates.clear();
    ghost_contact_pair_candidates.clear();
  }

  inline void
  clear_total_neighbor_list()
  {
    total_neighbor_list.clear();
  }


  // DEM containers
  typename dem_data_structures<dim>::cells_total_neighbor_list
    total_neighbor_list;

  // Modified by find_cell_neighbors class and shows the local neighbor cells of
  // all local cells in the triangulation
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;

  // Modified by find_cell_neighbors class and shows the ghost neighbor cells of
  // all local cells in the triangulation
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;

  // Map of vectors which contains all the local-local particle pairs in
  // adjacent cells which are collision candidates
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates;

  // Map of vectors which contains all the local-ghost particle pairs in
  // adjacent cells which are collision candidates
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_candidates;


  typename dem_data_structures<dim>::particle_floating_mesh_candidates
    particle_floating_mesh_candidates;

  // Container that contains all the contact information of particle-floating
  // mesh
  typename dem_data_structures<dim>::particle_floating_mesh_in_contact
    particle_floating_mesh_in_contact;
  typename dem_data_structures<dim>::particle_floating_wall_candidates
    particle_floating_wall_candidates;

  // Container that contains all the contact information of particle-floating
  // wall
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_floating_wall_in_contact;
  typename dem_data_structures<dim>::particle_wall_candidates
    particle_wall_candidates;

  // Container that contains all the contact information of particle-wall
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_in_contact;

  typename dem_data_structures<dim>::floating_mesh_information
    floating_mesh_info;

  typename dem_data_structures<dim>::particle_point_candidates
    particle_point_candidates;
  typename dem_data_structures<dim>::particle_line_candidates
    particle_line_candidates;

  // Container that contains all the contact information of particle-point
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_points_in_contact;

  // Container that contains all the contact information of particle-line
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_lines_in_contact;

  // Container that contains all the contact information of adjacent local-local
  // for calculation of the contact force of local-local particle pairs
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles;

  // Container that contains all the contact information of adjacent ghost-local
  // for calculation of the contact force of ghost-local particle pairs
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_adjacent_particles;

  // Container that contains the iterators to all local and ghost particles
  typename dem_data_structures<dim>::particle_index_iterator_map
    particle_container;

  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;

private:
  // Broad search objects
  ParticleParticleBroadSearch<dim>  particle_particle_broad_search_object;
  ParticleWallBroadSearch<dim>      particle_wall_broad_search_object;
  ParticlePointLineBroadSearch<dim> particle_point_line_broad_search_object;

  // Fine search objects
  ParticleParticleFineSearch<dim>  particle_particle_fine_search_object;
  ParticleWallFineSearch<dim>      particle_wall_fine_search_object;
  ParticlePointLineFineSearch<dim> particle_point_line_fine_search_object;

  // Other
  FindCellNeighbors<dim> cell_neighbors_object;
};



#endif // lethe_dem_container_manager_h
