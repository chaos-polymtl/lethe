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

#ifndef lethe_dem_container_manager_h
#  define lethe_dem_container_manager_h

template <int dim>
class DEMContainerManager
{
public:
  /**
   * Finds the neighbors lists of all the active cells in the input
   * triangulation.
   *
   * @param triangulation The triangulation to access the information of the cells
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   * @param has_floating_mesh A boolean to allow dealing with floating mesh neighbors
   *
   */

  void
  execute_cell_neighbors_search(
    const parallel::distributed::Triangulation<dim> &triangulation,
    const bool has_periodic_boundaries = false,
    const bool has_floating_mesh       = false);

  /**
   * Update the neighbors lists of all the active cells in the input
   * triangulation while clearing containers before the update.
   *
   * @param triangulation The triangulation to access the information of the cells
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   * @param has_floating_mesh A boolean to allow dealing with floating mesh neighbors
   *
   */

  void
  update_cell_neighbors(
    const parallel::distributed::Triangulation<dim> &triangulation,
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
                                  has_periodic_boundaries,
                                  has_floating_mesh);
  }

  /**
   * Manages to call update_fine_search_candidates() with contact data
   * containers to remove contact repetitions and to add new contact pairs to
   * the contact containers when particles are exchanged between processors.
   * Repeats the update for local particle-particle contacts, ghost
   * particle-particle, particle-wall contacts, particle-floating wall contacts
   * and particle-floating mesh contacts.
   *
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   */

  void
  update_contacts(const bool has_periodic_boundaries = false);

  /**
   * Updates the iterators to particles in local-local contact containers. This
   * is essential since sort_particles_into_subdomains_and_cells() and
   * exchange_ghost_particles() functions change the iterator to particles
   * everytime they are called.
   *
   * @tparam dim Dimensionality of the geometry which contains the particles
   * @param particle_handler The particle handler of particles
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   *
   */

  void
  update_local_particles_in_cells(
    const Particles::ParticleHandler<dim> &particle_handler,
    const bool                             has_periodic_boundaries = false);

  /**
   * Particle-particle broad search which finds the particle-particle contact
   * pairs candidates, local or ghost, which shows the collision pairs.
   * These collision pairs will be used in the fine search
   * to investigate if they are in contact or not.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   */

  void
  execute_particle_particle_broad_search(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    const bool                               has_periodic_boundaries = false);

  /**
   * Carries out the broad contact detection search using the
   * background triangulation for particle-walls contact.
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param boundary_cells_information Information of the boundary cells and
   * faces. This is the output of the FindBoundaryCellsInformation class
   * @param floating_wall Properties of the floating walls specified in the parameter handler file
   * @param simulation_time Simulation time
   * @param has_floating_mesh A boolean to allow dealing with floating mesh neighbors
   */

  void
  execute_particle_wall_broad_search(
    const Particles::ParticleHandler<dim> &           particle_handler,
    BoundaryCellsInformation<dim> &                   boundary_cell_object,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const bool has_floating_mesh = false);

  /**
   * Iterates over a vector of maps to see if the particles
   * which were in contact in the last time step, are still in contact or not.
   *
   * @param neighborhood_threshold A value which defines the neighbor particles
   * @param has_periodic_boundaries A boolean to allow periodic container manipulations
   * @param periodic_offset An offset used for tuning particle locations in periodic boundaries
   */

  void
  execute_particle_particle_fine_search(
    const double         neighborhood_threshold,
    const bool           has_periodic_boundaries = false,
    const Tensor<1, dim> periodic_offset         = Tensor<1, dim>());

  /**
   * Carries out the fine particle-wall contact detection
   *
   * @param floating_wall Properties of the floating walls specified in the parameter handler file
   * @param simulation_time Simulation time
   * @param neighborhood_threshold A value which defines the neighbor particles
   * @param has_floating_mesh A boolean to allow dealing with floating mesh neighbors
   */

  void
  execute_particle_wall_fine_search(
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_walls,
    const double                                      simulation_time,
    const double                                      neighborhood_threshold,
    const bool has_floating_mesh = false);

  /**
   * Carries out the storing of the floating mesh informations for background
   * and solid cells
   *
   * @param triangulation The triangulation to access the information of the cells
   * @param solids An object which manages a serial triangulation that represents a solid object.
   */

  void
  store_floating_mesh_info(
    const parallel::distributed::Triangulation<dim> &       triangulation,
    std::vector<std::shared_ptr<SerialSolid<dim - 1, dim>>> solids);


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
  typename dem_data_structures<dim>::floating_mesh_information
    floating_mesh_info;
  typename dem_data_structures<dim>::boundary_points_and_normal_vectors
    updated_boundary_points_and_normal_vectors;
  typename dem_data_structures<dim>::vector_on_boundary
    forces_boundary_information;
  typename dem_data_structures<dim>::vector_on_boundary
    torques_boundary_information;
  typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
    periodic_boundaries_cells_information;
  typename DEM::dem_data_structures<dim>::cell_container
    periodic_cells_container;

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

#endif // lethe_dem_container_manager_h
