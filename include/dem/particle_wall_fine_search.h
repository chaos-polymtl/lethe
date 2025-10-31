// SPDX-FileCopyrightText: Copyright (c) 2020, 2022, 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_fine_search_h
#define lethe_particle_wall_fine_search_h

#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

/**
 * @brief Iterate over the contact candidates from particle-wall broad search
 * (particle_wall_contact_pair_candidates) to add new contact pairs to the
 * particle_wall_pairs_in_contact container
 *
 * @param particle_wall_contact_pair_candidates The output of particle-wall broad search
 * which shows contact pair candidates
 * @param particle_wall_pairs_in_contact An unordered_map of maps which stores
 * all the particle-wall pairs which are physically in contact, and the
 * contact information in a struct. Note that the size of this unordered map
 * is equal to the number of particles
 */
template <int dim>
void
particle_wall_fine_search(
  const typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_pair_candidates,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_wall_pairs_in_contact);

/**
 * @brief Iterate over the contact candidates from particle-floating wall broad
 * search (particle_floating_wall_candidates) to add new contact pairs to the
 * particle_floating_wall_pairs_in_contact container
 *
 * @param particle_floating_wall_contact_candidates The output of particle-floating wall
 * broad search which shows contact pair candidates
 * @param floating_wall_properties Properties of floating walls defined in the
 * parameter handler
 * @param simulation_time Simulation time
 * @param particle_floating_wall_in_contact An unordered_map of maps which stores
 * all the particle-floating wall pairs which are in contact, and
 * the contact information in a struct. Note that the size of this unordered
 * map is equal to the number of particles
 */
template <int dim>
void
particle_floating_wall_fine_search(
  const typename DEM::dem_data_structures<
    dim>::particle_floating_wall_candidates
    &particle_floating_wall_contact_candidates,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_wall_in_contact
    &particle_floating_wall_in_contact);

/**
 * @brief Iterate over the contact candidates from particle-floating mesh broad
 * search (particle_floating_mesh_contact_candidates) to add new contact pairs
 * to the particle_floating_mesh_potentially_in_contact container
 *
 * @param particle_floating_mesh_contact_candidates The output of particle-floating mesh
 * broad search which shows contact pair candidates
 * @param particle_floating_mesh_potentially_in_contact A map of maps which stores
 * all the particle-floating mesh pairs which are in contact
 */
template <int dim>
void
particle_floating_mesh_fine_search(
  const typename DEM::dem_data_structures<
    dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<
    dim>::particle_floating_mesh_potentially_in_contact
    &particle_floating_mesh_potentially_in_contact);

#endif
