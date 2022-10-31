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

#include <dem/data_containers.h>

using namespace DEM;

#ifndef lethe_dem_container_manager_h
#  define lethe_dem_container_manager_h

template <int dim>
class DEMContainerManager
{
public:
  // DEM containers
  typename dem_data_structures<dim>::cells_total_neighbor_list
    total_neighbor_list;
  typename dem_data_structures<dim>::floating_mesh_information
    floating_mesh_info;

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
};



#endif // lethe_dem_container_manager_h
