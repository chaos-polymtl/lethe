// SPDX-FileCopyrightText: Copyright (c) 2020-2022, 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_wall_broad_search_h
#define lethe_particle_wall_broad_search_h

#include <dem/adaptive_sparse_contacts.h>
#include <dem/boundary_cells_info_struct.h>
#include <dem/contact_info.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/particles/particle_handler.h>

#include <map>

using namespace dealii;

/**
 * @brief Finds unordered map of tuples (tuple of particle located in
 * boundary cells, normal vector of the boundary face, a
 * point on the face and the corresponding boundary cell) which shows the
 * candidate particle-wall collision pairs. These collision candidates will be
 * investigated in the fine search to check if they are in contact or not.
 *
 * @param boundary_cells_information Information of the boundary cells and
 * faces. This is the output of the FindBoundaryCellsInformation class.
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param particle_wall_contact_candidates A two-layered unordered map of tuples. Each
 * tuple contains a particle located near boundaries, the normal vector of
 * the corresponding face boundary, a point on the boundary and the boundary
 * cell. The contact pair is used in the fine search.
 */
template <int dim>
void
find_particle_wall_contact_pairs(
  const std::map<int, boundary_cells_info_struct<dim>>
                                        &boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_candidates);

/**
 * @brief Finds unordered map of tuples (tuple of particle located in
 * boundary cells, normal vector of the boundary face, a
 * point on the face and the corresponding boundary cell) which shows the
 * candidate particle-wall collision pairs. These collision candidates will be
 * investigated in the fine search to check if they are in contact or not.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param boundary_cells_information Information of the boundary cells and
 * faces. This is the output of the FindBoundaryCellsInformation class.
 * @param particle_handler Particle handler of particles located in boundary
 * cells.
 * @param particle_wall_contact_candidates A two-layered unordered map of tuples. Each
 * tuple contains a particle located near boundaries, the normal vector of
 * the corresponding face boundary, a point on the boundary and the boundary
 * cell. The contact pair is used in the fine search.
 * @param sparse_contacts_object The object that contains the
 * information about the mobility status of cells
 */
template <int dim, typename PropertiesIndex>
void
find_particle_wall_contact_pairs(
  const std::map<int, boundary_cells_info_struct<dim>>
                                        &boundary_cells_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &particle_wall_contact_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

/**
 * @brief Find a two-layered unordered map of particle iterators which shows the
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
 * @param particle_floating_wall_candidates Output of particle-floating wall
 * broad search which contains all the particle-floating wall collision
 * candidates
 */
template <int dim>
void
find_particle_floating_wall_contact_pairs(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
                                        &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &particle_floating_wall_candidates);

/**
 * @brief Find a two-layered unordered map of particle iterators which shows the
 * candidate particle-floating wall collision candidates. These collision
 * pairs will be investigated in the fine search to check if they are in
 * contact or not
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param boundary_cells_for_floating_walls Boundary cells located adjacent to
 * floating walls
 * @param particle_handler Particle handler of particles located in boundary
 * cells
 * @param floating_wall_properties Properties of the floating walls specified
 * in the parameter handler file
 * @param simulation_time Simulation time
 * @param particle_floating_wall_candidates Output of particle-floating wall
 * broad search which contains all the particle-floating wall collision
 * candidates
 * @param sparse_contacts_object The object that contains the
 * information about the mobility status of cells
 */
template <int dim, typename PropertiesIndex>
void
find_particle_floating_wall_contact_pairs(
  const std::unordered_map<
    types::global_dof_index,
    std::set<typename Triangulation<dim>::active_cell_iterator>>
                                        &boundary_cells_for_floating_walls,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
  const double                                      simulation_time,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &particle_floating_wall_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

/**
 * @brief Find a two-layered unordered map
 * (particle_floating_mesh_contact_candidates) of particle iterators that
 * shows the candidate particle-floating mesh collision candidates. These
 * collision pairs will be investigated in the fine search to check if they
 * are in contact or not
 *
 * @param solid_surfaces_mesh_information Information of the solid surfaces mapped
 * in the background triangulation.
 * @param particle_handler
 * @param particle_floating_mesh_contact_candidates Particle-floating mesh contact
 * candidates.
 * @param cells_total_neighbor_list A container in which all the neighbor cells
 * of the local cells are stored.
 */
template <int dim>
void
particle_solid_surfaces_contact_search(
  const typename DEM::dem_data_structures<dim>::solid_surfaces_mesh_information
                                        &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
    &cells_total_neighbor_list);


/**
 * @brief Find a two-layered unordered map
 * (particle_floating_mesh_contact_candidates) of particle iterators that
 * shows the candidate particle-floating mesh collision candidates. These
 * collision pairs will be investigated in the fine search to check if they
 * are in contact or not
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param solid_surfaces_mesh_information Information of the solid surfaces mapped
 * in the background triangulation.
 * @param particle_handler
 * @param particle_floating_mesh_contact_candidates Particle-floating mesh contact
 * candidates
 * @param cells_total_neighbor_list A container in which all the neighbor cells
 * of the local cells are stored
 * @param sparse_contacts_object The object that contains the
 * information about the mobility status of cells
 */
template <int dim, typename PropertiesIndex>
void
particle_solid_surfaces_contact_search(
  const typename DEM::dem_data_structures<dim>::solid_surfaces_mesh_information
                                        &solid_surfaces_mesh_information,
  const Particles::ParticleHandler<dim> &particle_handler,
  typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
    &particle_floating_mesh_contact_candidates,
  typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
                                                     &cells_total_neighbor_list,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);


/**
 * @brief Store the candidate particle-wall contact pairs.
 *
 * @param particle_iterator The particle iterator to start boundary id.
 * @param boundary_cells_content The info of the current boundary cell to
 * store.
 * @param contact_pair_candidates A map which will contain all the
 * particle/wall pairs candidate.
 */
template <int dim>
inline void
store_candidates(
  const typename Particles::ParticleHandler<
    dim>::particle_iterator_range::iterator &particle_iterator,
  const boundary_cells_info_struct<dim>     &boundary_cells_content,
  typename DEM::dem_data_structures<dim>::particle_wall_candidates
    &contact_pair_candidates)
{
  // Find the contact candidate container of the particle
  const types::particle_index particle_id = particle_iterator->get_id();
  auto candidates_container_it = contact_pair_candidates.find(particle_id);

  // Reserve arbitrary vector capacity and store if the particle does not have
  // contact candidate yet
  if (candidates_container_it == contact_pair_candidates.end())
    {
      auto pair_it_bool = contact_pair_candidates.emplace(
        particle_id,
        ankerl::unordered_dense::map<
          DEM::global_face_id,
          std::tuple<Particles::ParticleIterator<dim>,
                     Tensor<1, dim>,
                     Point<dim>,
                     DEM::global_face_id>>());

      candidates_container_it = pair_it_bool.first;
    }

  // Store particle ids from the selected particle iterator
  candidates_container_it->second.emplace(
    boundary_cells_content.global_face_id,
    std::make_tuple(particle_iterator,
                    boundary_cells_content.normal_vector,
                    boundary_cells_content.point_on_face,
                    boundary_cells_content.boundary_id));
}

/**
 * @brief Store the candidate particle-floating wall contact pairs.
 *
 * @param particle_iterator The particle iterator to start boundary id.
 * @param floating_wall_id The floating wall id.
 * @param contact_pair_candidates A map which will contain all the
 * particle/wall pairs candidate.
 */
template <int dim>
inline void
store_candidates(
  const typename Particles::ParticleHandler<
    dim>::particle_iterator_range::iterator &particle_iterator,
  const DEM::global_face_id                 &floating_wall_id,
  typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
    &contact_pair_candidates);

#endif
