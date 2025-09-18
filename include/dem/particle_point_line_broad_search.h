// SPDX-FileCopyrightText: Copyright (c) 2020, 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_point_line_broad_search_h
#define lethe_particle_point_line_broad_search_h

#include <dem/adaptive_sparse_contacts.h>
#include <dem/contact_info.h>
#include <dem/data_containers.h>

#include <deal.II/particles/particle_iterator.h>

using namespace dealii;

/**
 * @brief Find a map of pairs (pair of particle and the boundary vertex
 * location) which shows the candidate particle-point collision pairs. These
 * collision pairs will be investigated in the fine search to check if they
 * are in contact or not.
 *
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param boundary_cells_with_points A container of cells which are located at
 * boundaries with only one vertex.
 * @param particle_point_contact_candidates The map of particle-point pairs.
 * Each element of map (pair) contains a contact pair particle located near
 * boundaries with vertices and the vertex location.
 */
template <int dim>
void
find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_point_info<dim>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<dim>::particle_point_candidates
    &particle_point_contact_candidates);

/**
 * @brief Find a map of pairs (pair of particle and the boundary vertex
 * location) which shows the candidate particle-point collision pairs. These
 * collision pairs will be investigated in the fine search to check if they
 * are in contact or not.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param boundary_cells_with_points A container of cells which are located at
 * boundaries with only one vertex.
 * @param particle_point_contact_candidates The map of particle-point pairs.
 * Each element of map (pair) contains a contact pair particle located near
 * boundaries with vertices and the vertex location.
 * @param sparse_contacts_object The Adaptive Sparse Contacts for mobility
 * status checks.
 */
template <int dim, typename PropertiesIndex>
void
find_particle_point_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_point_info<dim>>
    &boundary_cells_with_points,
  typename DEM::dem_data_structures<dim>::particle_point_candidates
    &particle_point_contact_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

/**
 * @brief Find a map of tuples (tuple of particle and the locations of
 * beginning and ending vertices of the boundary lines) which shows the
 * candidate particle-line collision pairs. These collision pairs will be
 * investigated in the fine search to check if they are in contact or not.
 *
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param boundary_cells_with_lines A container of cells which are located at
 * boundaries with only one line.
 * @param particle_line_contact_candidates Each element of map (tuple) contains
 * a particle and the locations of beginning and ending vertices of the boundary
 * lines.
 */
template <int dim>
void
find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<dim>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<dim>::particle_line_candidates
    &particle_line_contact_candidates);

/**
 * @brief Find a map of tuples (tuple of particle and the locations of
 * beginning and ending vertices of the boundary lines) which shows the
 * candidate particle-line collision pairs. These collision pairs will be
 * investigated in the fine search to check if they are in contact or not.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param boundary_cells_with_lines A container of cells which are located at
 * boundaries with only one line.
 * @param particle_line_contact_candidates Each element of map (tuple) contains
 * a particle and the locations of beginning and ending vertices of the boundary
 * lines.
 * @param sparse_contacts_object The Adaptive Sparse Contacts for mobility
 * status checks.
 */
template <int dim, typename PropertiesIndex>
void
find_particle_line_contact_pairs(
  const Particles::ParticleHandler<dim> &particle_handler,
  const std::unordered_map<std::string, cell_line_info<dim>>
    &boundary_cells_with_lines,
  typename DEM::dem_data_structures<dim>::particle_line_candidates
    &particle_line_contact_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

#endif
