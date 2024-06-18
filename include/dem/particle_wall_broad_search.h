/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2024 by the Lethe authors
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

#include <dem/adaptive_sparse_contacts.h>
#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_solver_parameters.h>

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
 * Broad search is used to obtain all the particles located at boundary cells
 */

template <int dim>
class ParticleWallBroadSearch
{
public:
  ParticleWallBroadSearch();

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
                                          &boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates
      &particle_wall_contact_candidates);

  void
  find_particle_wall_contact_pairs(
    const std::map<int, boundary_cells_info_struct<dim>>
                                          &boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_wall_candidates
                                      &particle_wall_contact_candidates,
    const AdaptiveSparseContacts<dim> &sparse_contacts_object);

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
      types::global_dof_index,
      std::set<typename Triangulation<dim>::active_cell_iterator>>
                                          &boundary_cells_for_floating_walls,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double                                      simulation_time,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
      &particle_floating_wall_candidates);

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
    const AdaptiveSparseContacts<dim> &sparse_contacts_object);

  /**
   * Finds a two-layered unordered map
   * (particle_floating_mesh_contact_candidates) of particle iterators that
   * shows the candidate particle-floating mesh collision candidates. These
   * collision pairs will be investigated in the fine search to check if they
   * are in contact or not
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

  void
  particle_solid_surfaces_contact_search(
    const typename DEM::dem_data_structures<
      dim>::solid_surfaces_mesh_information &solid_surfaces_mesh_information,
    const Particles::ParticleHandler<dim>   &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
      &particle_floating_mesh_contact_candidates,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
      &cells_total_neighbor_list);

  void
  particle_solid_surfaces_contact_search(
    const typename DEM::dem_data_structures<
      dim>::solid_surfaces_mesh_information &solid_surfaces_mesh_information,
    const Particles::ParticleHandler<dim>   &particle_handler,
    typename DEM::dem_data_structures<dim>::particle_floating_mesh_candidates
      &particle_floating_mesh_contact_candidates,
    typename DEM::dem_data_structures<dim>::cells_total_neighbor_list
                                      &cells_total_neighbor_list,
    const AdaptiveSparseContacts<dim> &sparse_contacts_object);

private:
  /**
   * @brief Stores the candidate particle-wall contact pairs.
   *
   * @param particle_iterator The particle iterator to start boundary id.
   * @param boundary_cells_content The info of the current boundary cell to
   * store.
   * @param contact_pair_candidates A map which will contain all the
   * particle/wall pairs candidate.
   */
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
   * @brief Stores the candidate particle-floating wall contact pairs.
   *
   * @param particle_iterator The particle iterator to start boundary id.
   * @param floating_wall_id The floating wall id.
   * @param contact_pair_candidates A map which will contain all the
   * particle/wall pairs candidate.
   */
  inline void
  store_candidates(
    const typename Particles::ParticleHandler<
      dim>::particle_iterator_range::iterator &particle_iterator,
    const DEM::global_face_id                 &floating_wall_id,
    typename DEM::dem_data_structures<dim>::particle_floating_wall_candidates
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
          std::unordered_map<DEM::global_face_id,
                             Particles::ParticleIterator<dim>>());

        candidates_container_it = pair_it_bool.first;
      }

    // Store particle ids from the selected particle iterator
    candidates_container_it->second.emplace(floating_wall_id,
                                            particle_iterator);
  }
};

#endif /* particle_wall_broad_search_h */
