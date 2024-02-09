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
#include <dem/disable_contacts.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/find_cell_neighbors.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_particle_fine_search.h>
#include <dem/particle_point_line_broad_search.h>
#include <dem/particle_point_line_fine_search.h>
#include <dem/particle_wall_broad_search.h>
#include <dem/particle_wall_fine_search.h>
#include <dem/update_fine_search_candidates.h>
#include <dem/update_local_particle_containers.h>

#include <deal.II/particles/particle_handler.h>

using namespace DEM;

#ifndef lethe_dem_contact_manager_h
#  define lethe_dem_contact_manager_h

/**
 * @brief Manage the numerous of contact detection search operations in the DEM.
 *
 * This class mostly calls proper functions in regards the type of contacts for
 * the contact detection and updates of data containers.
 */
template <int dim>
class DEMContactManager
{
public:
  /**
   * @brief Execute functions to find the lists of cell neighbors of the active
   * cells in the triangulation.
   *
   * For contacts of particles, the neighbor search excludes the repetitions of
   * neighbors: B is neighbor of A but A will not be neighbor of B. For contacts
   * of particles with floating meshes, the neighbor search includes all
   * neighbors. These floating mesh contacts searches need all the particles
   * located in all the neighbor cells of the main cell to detect
   * collisions with the floating mesh.
   *
   * @param[in] triangulation The triangulation to access the information of the
   * cells
   * @param[in] periodic_boundaries_cells_information Information of periodic
   * cells used if periodic boundaries are enabled (next parameter).
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   * @param[in] has_floating_mesh Allow the computation of full neighbor lists
   * of active cells computation if required.
   */
  void
  execute_cell_neighbors_search(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
               periodic_boundaries_cells_information,
    const bool has_periodic_boundaries = false,
    const bool has_floating_mesh       = false);

  /**
   * @brief Execute functions to clear and update the neighbor lists of all the
   * active cells in the triangulation.
   *
   * @param[in] triangulation The triangulation to access the information of the
   * cells
   * @param[in] periodic_boundaries_cells_information Information of periodic
   * cells used if periodic boundaries are enabled (next parameter).
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   * @param[in] has_floating_mesh Allow the computation of full neighbor lists
   * of active cells computation if required.
   */
  void
  update_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
               periodic_boundaries_cells_information,
    const bool has_periodic_boundaries = false,
    const bool has_floating_mesh       = false)
  {
    cells_local_neighbor_list.clear();
    cells_ghost_neighbor_list.clear();
    cells_local_periodic_neighbor_list.clear();
    cells_ghost_periodic_neighbor_list.clear();
    cells_ghost_local_periodic_neighbor_list.clear();
    total_neighbor_list.clear();

    execute_cell_neighbors_search(triangulation,
                                  periodic_boundaries_cells_information,
                                  has_periodic_boundaries,
                                  has_floating_mesh);
  }

  /**
   * @brief Execute functions to update the contact data containers with the new
   * contact pairs.
   *
   * Call proper functions to remove contact repetitions and to add new contact
   * pairs to the contact containers when particles are exchanged between
   * processors after the fine search.
   * Contact pairs are local particle-particle, ghost particle-particle,
   * particle-wall, particle-floating wall contacts and particle-floating mesh
   * contacts.
   *
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   */
  void
  update_contacts(const bool has_periodic_boundaries = false);

  /**
   * @breif Execute functions to update the particle iterators in local-local
   * contact containers.
   *
   * This is essential since sort_particles_into_subdomains_and_cells() and
   * exchange_ghost_particles() functions change the particle iterators
   * everytime they are called.
   *
   * @param[in] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] clear_contact_structures Allow clearing the contact structures.
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   */
  void
  update_local_particles_in_cells(
    const Particles::ParticleHandler<dim> &particle_handler,
    const bool                             clear_contact_structures,
    const bool                             has_periodic_boundaries = false);

  /**
   * @brief Execute the particle-particle broad searches.
   *
   * Calls proper functions to find the candidates of local and ghost
   * particle-particle contact pairs and the periodic particle-particle contacts
   * if required. These contact pairs will be used in the fine search step to
   * investigate if they are in contact.
   *
   * @param[in] particle_handler Storage of particles and their accessor
   * functions. search[in]
   * @param has_periodic_boundaries Allow manipulations of periodic containers
   * if required.
   */
  void
  execute_particle_particle_broad_search(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const bool                               has_periodic_boundaries = false);

  /**
   * @brief Execute the particle-particle broad searches with disabling of
   * contacts.
   *
   * Calls proper functions to find the candidates of local and ghost
   * particle-particle contact pairs and the periodic particle-particle contacts
   * if required. These contact pairs will be used in the fine search step to
   * investigate if they are in contact. This version of the function is used
   * when disabling particle contacts regards mobility is enabled.
   *
   * @param[in,out] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] disable_particle_contact_object Allow to check the mobility
   * status of cells
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   */
  void
  execute_particle_particle_broad_search(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const DisableContacts<dim>              &disable_particle_contact_object,
    const bool                               has_periodic_boundaries = false);

  /**
   * @brief Execute the particle-wall broad searches.
   *
   * Calls proper functions to find the candidates of particle-wall contact
   * pairs and the particle-floating wall or particle-floating mesh if required.
   * These contact pairs will be used in the fine search step to investigate if
   * they are in contact.
   *
   * @param[in] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] boundary_cells_object Information of the boundary cells and
   * faces.
   * @param[in] floating_mesh_info Mapping of floating meshes.
   * @param[in] floating_wall Properties of the floating walls.
   * @param[in] simulation_time Current simulation time.
   * @param[in] has_floating_mesh Allow dealing with floating mesh neighbors.
   */
  void
  execute_particle_wall_broad_search(
    const Particles::ParticleHandler<dim> &particle_handler,
    BoundaryCellsInformation<dim>         &boundary_cell_object,
    const typename dem_data_structures<dim>::floating_mesh_information
                                                      floating_mesh_info,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const bool has_floating_mesh = false);

  /**
   * @brief Execute the particle-wall broad searches with disabling of contacts.
   *
   * Calls proper functions to find the candidates of particle-wall contact
   * pairs and the particle-floating wall or particle-floating mesh if required.
   * These contact pairs will be used in the fine search step to investigate if
   * they are in contact. This version of the function is used when disabling
   * particle contacts regards mobility is enabled.
   *
   * @param[in] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] boundary_cells_object Information of the boundary cells and
   * faces.
   * @param[in] floating_mesh_info Mapping of floating meshes.
   * @param[in] floating_wall Properties of the floating walls.
   * @param[in] simulation_time Current simulation time.
   * @param[in] disable_particle_contact_object Allow to check the mobility
   * status of cells
   * @param[in] has_floating_mesh Allow dealing with floating mesh neighbors.
   */
  void
  execute_particle_wall_broad_search(
    const Particles::ParticleHandler<dim> &particle_handler,
    BoundaryCellsInformation<dim>         &boundary_cell_object,
    const typename dem_data_structures<dim>::floating_mesh_information
                                                      floating_mesh_info,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const DisableContacts<dim> &disable_particle_contact_object,
    const bool                  has_floating_mesh = false);

  /**
   * @brief Execute the particle-particles fine searches.
   *
   * Executes functions that update the particle contacts pairs containers and
   * compute the contact information of the collision pairs.
   *
   * @param[in] neighborhood_threshold Threshold value of contact detection.
   * @param[in] has_periodic_boundaries Allow manipulations of periodic
   * containers if required.
   * @param[in] periodic_offset Offset used for tuning particle locations in
   * periodic boundaries.
   */
  void
  execute_particle_particle_fine_search(
    const double         neighborhood_threshold,
    const bool           has_periodic_boundaries = false,
    const Tensor<1, dim> periodic_offset         = Tensor<1, dim>());

  /**
   * @brief Execute the particle-wall fine searches.
   *
   * Executes functions that update the particle-wall contacts pairs containers
   * and compute the contact information of the collision particle-wall.
   *
   * @param[in] floating_wall Properties of the floating walls.
   * @param[in] simulation_time Current simulation time.
   * @param[in] neighborhood_threshold Threshold value of contact detection.
   * @param[in] has_floating_mesh Allow the fine search with floating meshes.
   */
  void
  execute_particle_wall_fine_search(
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const double                                      neighborhood_threshold,
    const bool has_floating_mesh = false);


  // Container with the iterators to all local and ghost particles
  typename dem_data_structures<dim>::particle_index_iterator_map
    particle_container;

  // Container that shows the local/ghost neighbor cells of all local cells in
  // the triangulation
  typename dem_data_structures<dim>::cells_total_neighbor_list
    total_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_local_periodic_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_periodic_neighbor_list;
  typename dem_data_structures<dim>::cells_neighbor_list
    cells_ghost_local_periodic_neighbor_list;

  // Container with all collision candidate particles within adjacent cells
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_local_contact_pair_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    local_contact_pair_periodic_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_contact_pair_periodic_candidates;
  typename dem_data_structures<dim>::particle_particle_candidates
    ghost_local_contact_pair_periodic_candidates;
  typename dem_data_structures<dim>::particle_floating_mesh_candidates
    particle_floating_mesh_candidates;
  typename dem_data_structures<dim>::particle_floating_wall_candidates
    particle_floating_wall_candidates;
  typename dem_data_structures<dim>::particle_point_candidates
    particle_point_candidates;
  typename dem_data_structures<dim>::particle_line_candidates
    particle_line_candidates;
  typename dem_data_structures<dim>::particle_wall_candidates
    particle_wall_candidates;

  // Container with all the contact information of the object in contact
  // with a particle
  typename dem_data_structures<dim>::particle_floating_mesh_in_contact
    particle_floating_mesh_in_contact;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_floating_wall_in_contact;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_in_contact;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_points_in_contact;
  typename dem_data_structures<dim>::particle_point_line_contact_info
    particle_lines_in_contact;

  // Container with all the contact information of adjacent
  // local/ghost-local for pairwise contact force calculation
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_periodic_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_periodic_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_local_periodic_adjacent_particles;

  // Containers with other information
  typename DEM::dem_data_structures<dim>::cell_vector periodic_cells_container;

private:
  // Broad search objects
  ParticleParticleBroadSearch<dim>  particle_particle_broad_search_object;
  ParticleWallBroadSearch<dim>      particle_wall_broad_search_object;
  ParticlePointLineBroadSearch<dim> particle_point_line_broad_search_object;

  // Fine search objects
  ParticleParticleFineSearch<dim>  particle_particle_fine_search_object;
  ParticleWallFineSearch<dim>      particle_wall_fine_search_object;
  ParticlePointLineFineSearch<dim> particle_point_line_fine_search_object;

  // Other relevant objects
  FindCellNeighbors<dim> cell_neighbors_object;
};

#endif // lethe_dem_contact_manager_h
