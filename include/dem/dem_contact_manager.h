// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_dem_contact_manager_h
#define lethe_dem_contact_manager_h


#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_particle_broad_search.h>
#include <dem/particle_wall_broad_search.h>


using namespace DEM;

/**
 * @brief Manage the numerous of contact detection search operations in the DEM.
 *
 * This class mostly calls proper functions in regards the type of contacts for
 * the contact detection and updates of data containers.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
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
   * cells.
   * @param[in] periodic_boundaries_cells_information Information of periodic
   * cells used if periodic boundaries are enabled (next parameter).
   */
  void
  execute_cell_neighbors_search(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      periodic_boundaries_cells_information);

  /**
   * @brief Execute functions to clear and update the neighbor lists of all the
   * active cells in the triangulation.
   *
   * @param[in] triangulation The triangulation to access the information of the
   * cells.
   * @param[in] periodic_boundaries_cells_information Information of periodic
   * cells used if periodic boundaries are enabled (next parameter).
   */
  void
  update_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      periodic_boundaries_cells_information)
  {
    cells_local_neighbor_list.clear();
    cells_ghost_neighbor_list.clear();
    cells_local_periodic_neighbor_list.clear();
    cells_ghost_periodic_neighbor_list.clear();
    cells_ghost_local_periodic_neighbor_list.clear();
    total_neighbor_list.clear();

    execute_cell_neighbors_search(triangulation,
                                  periodic_boundaries_cells_information);
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
   */
  void
  update_contacts();

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
   */
  void
  update_local_particles_in_cells(
    const Particles::ParticleHandler<dim> &particle_handler);

  /**
   * @brief Execute the particle-particle broad searches with adaptive sparse
   * contacts.
   *
   * Calls proper functions to find the candidates of local and ghost
   * particle-particle contact pairs and the periodic particle-particle contacts
   * if required. These contact pairs will be used in the fine search step to
   * investigate if they are in contact.
   * It checks if the adaptive sparse contacts is enabled and use proper
   * functions.
   *
   * @param[in,out] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] sparse_particle_contact_object Allow to check the mobility
   * status of cells.
   */
  void
  execute_particle_particle_broad_search(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const AdaptiveSparseContacts<dim, PropertiesIndex>
      &sparse_particle_contact_object);

  /**
   * @brief Execute the particle-wall broad searches with adaptive sparse
   * contacts.
   *
   * Calls proper functions to find the candidates of particle-wall contact
   * pairs and the particle-floating wall or particle-floating mesh if required.
   * These contact pairs will be used in the fine search step to investigate if
   * they are in contact.
   * It checks if the adaptive sparse contacts is enabled and use proper
   * functions.
   *
   * @param[in] particle_handler Storage of particles and their accessor
   * functions.
   * @param[in] boundary_cells_object Information of the boundary cells and
   * faces.
   * @param[in] solid_surfaces_mesh_info Mapping of solid surfaces meshes.
   * @param[in] floating_wall Properties of the floating walls.
   * @param[in] simulation_time Current simulation time.
   * @param[in] sparse_particle_contact_object Allow to check the mobility
   * status of cells.
   */
  void
  execute_particle_wall_broad_search(
    const Particles::ParticleHandler<dim> &particle_handler,
    BoundaryCellsInformation<dim>         &boundary_cell_object,
    const typename dem_data_structures<dim>::solid_surfaces_mesh_information
                                                      solid_surfaces_mesh_info,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const AdaptiveSparseContacts<dim, PropertiesIndex>
      &sparse_particle_contact_object);

  /**
   * @brief Execute the particle-particles fine searches.
   *
   * Executes functions that update the particle contacts pairs containers and
   * compute the contact information of the collision pairs.
   *
   * @param[in] neighborhood_threshold Threshold value of contact detection.
   */
  void
  execute_particle_particle_fine_search(const double neighborhood_threshold);

  /**
   * @brief Execute the particle-wall fine searches.
   *
   * Executes functions that update the particle-wall contacts pairs containers
   * and compute the contact information of the collision particle-wall.
   *
   * @param[in] floating_wall Properties of the floating walls.
   * @param[in] simulation_time Current simulation time.
   * @param[in] neighborhood_threshold Threshold value of contact detection.
   */
  void
  execute_particle_wall_fine_search(
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const double                                      neighborhood_threshold);

  /**
   * @brief Set the constant periodic offset for the periodic boundaries. If
   * they are no periodic boundaries, the default offset is zeros;
   *
   * @param[in] periodic_offset Offset used for tuning particle locations in
   * periodic boundaries.
   */
  inline void
  set_periodic_offset(const Tensor<1, dim> &offset)
  {
    this->periodic_offset = offset;
  }

  /**
   * @brief Return the particle-floating mesh contact container.
   */
  inline typename dem_data_structures<
    dim>::particle_floating_mesh_potentially_in_contact &
  get_particle_floating_mesh_potentially_in_contact()
  {
    return particle_floating_mesh_potentially_in_contact;
  }

  /**
   * @brief Return the particle-floating wall contact container.
   */
  inline typename dem_data_structures<dim>::particle_wall_in_contact &
  get_particle_floating_wall_in_contact()
  {
    return particle_floating_wall_in_contact;
  }

  /**
   * @brief Return the particle-wall contact container.
   */
  inline typename dem_data_structures<dim>::particle_wall_in_contact &
  get_particle_wall_in_contact()
  {
    return particle_wall_in_contact;
  }

  /**
   * @brief Return the particle-line contact container.
   */
  inline typename dem_data_structures<dim>::particle_line_candidates &
  get_particle_lines_in_contact()
  {
    return particle_lines_in_contact;
  }

  /**
   * @brief Return the particle-point contact container.
   */
  inline typename dem_data_structures<dim>::particle_point_candidates &
  get_particle_points_in_contact()
  {
    return particle_points_in_contact;
  }

  /**
   * @brief Return the local-local particle contact container.
   */
  inline typename dem_data_structures<dim>::adjacent_particle_pairs &
  get_local_adjacent_particles()
  {
    return local_adjacent_particles;
  }

  /**
   * @brief Return the local-ghost particle contact container.
   */
  inline typename dem_data_structures<dim>::adjacent_particle_pairs &
  get_ghost_adjacent_particles()
  {
    return ghost_adjacent_particles;
  }

  /**
   * @brief Return the local-local periodic particle contact container.
   */
  inline typename dem_data_structures<dim>::adjacent_particle_pairs &
  get_local_local_periodic_adjacent_particles()
  {
    return local_local_periodic_adjacent_particles;
  }

  /**
   * @brief Return the local-ghost periodic particle contact container.
   */
  inline typename dem_data_structures<dim>::adjacent_particle_pairs &
  get_local_ghost_periodic_adjacent_particles()
  {
    return local_ghost_periodic_adjacent_particles;
  }

  /**
   * @brief Return the ghost-local periodic particle contact container.
   */
  inline typename dem_data_structures<dim>::adjacent_particle_pairs &
  get_ghost_local_periodic_adjacent_particles()
  {
    return ghost_local_periodic_adjacent_particles;
  }

  /**
   * @brief Return the local particle-particle contact candidates.
   */
  inline typename dem_data_structures<dim>::particle_particle_candidates &
  get_local_contact_pair_candidates()
  {
    return local_contact_pair_candidates;
  }

  /**
   * @brief Return the ghost particle-particle contact candidates.
   */
  inline typename dem_data_structures<dim>::particle_particle_candidates &
  get_ghost_contact_pair_candidates()
  {
    return ghost_contact_pair_candidates;
  }

private:
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
  typename dem_data_structures<
    dim>::particle_floating_mesh_potentially_in_contact
    particle_floating_mesh_potentially_in_contact;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_floating_wall_in_contact;
  typename dem_data_structures<dim>::particle_wall_in_contact
    particle_wall_in_contact;
  typename dem_data_structures<dim>::particle_line_candidates
    particle_lines_in_contact;
  typename dem_data_structures<dim>::particle_point_candidates
    particle_points_in_contact;


  // Container with all the contact information of adjacent
  // local/ghost-local for pairwise contact force calculation
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_local_periodic_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    local_ghost_periodic_adjacent_particles;
  typename dem_data_structures<dim>::adjacent_particle_pairs
    ghost_local_periodic_adjacent_particles;

  // Containers with other information
  typename DEM::dem_data_structures<dim>::cell_vector periodic_cells_container;

private:
  Tensor<1, dim> periodic_offset = Tensor<1, dim>();
};

#endif
