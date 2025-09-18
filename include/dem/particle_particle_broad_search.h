// SPDX-FileCopyrightText: Copyright (c) 2020-2022, 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_particle_particle_broad_search_h
#define lethe_particle_particle_broad_search_h

#include <dem/adaptive_sparse_contacts.h>
#include <dem/data_containers.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

/**
 * @brief Finds a vector of pairs (particle_particle_candidates) which shows the
 * candidate particle-particle collision pairs. These collision pairs will be
 * used in the fine search to investigate if they are in contact or not.
 *
 * @param[in] particle_handler The particle handler of particles in the broad
 * search
 * @param[in] cells_local_neighbor_list  A vector (with size equal to the number
 * of local cells) of vectors. Each sub-vector have a size equal to the number
 * of adjacent local cells plus one. First element of each sub-vector shows the
 * main cell itself.
 * @param[in] cells_ghost_neighbor_list A vector (with size equal to the number
 * of local cells) of vectors. Each sub-vector have a size equal to the number
 * of adjacent ghost cells plus one. First element of each sub-vector shows the
 * main cell itself.
 * @param[out] local_contact_pair_candidates Ankerl unordered dense map. Stores
 * potential pairs of local-local particle in contact without redundancy.
 * Keys are particle ids and mapped types are vectors of particle ids.
 * @param[out] ghost_contact_pair_candidates Ankerl unordered dense map. Stores
 * potential pairs of local-ghost particle in contact. Keys are particle ids and
 * mapped types are vectors of particle ids.
 */
template <int dim>
void
find_particle_particle_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_candidates);

/**
 * @brief Finds a vector of pairs (particle_particle_candidates) which shows the
 * candidate particle-particle collision pairs. These collision pairs will be
 * used in the fine search to investigate if they are in contact or not.
 * This version of the function is used when adaptive sparse contacts is
 * enabled.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 *
 * @param particle_handler The particle handler of particles in the broad
 * search
 * @param[in] cells_local_neighbor_list  A vector (with size equal to the number
 * of local cells) of vectors. Each sub-vector have a size equal to the number
 * of adjacent local cells of the main cell plus one. The first element of each
 * sub-vector shows the main cell itself.
 * @param[in] cells_ghost_neighbor_list A vector (with size equal to the number
 * of local cells) of vectors. Each sub-vector have a size equal to the number
 * of adjacent ghost cells of the main cell plus one. The first element of each
 * sub-vector shows the main cell itself.
 * @param[out] local_contact_pair_candidates Ankerl unordered dense map. Stores
 * potential pairs of local-local particles in contact without redundancy.
 * Keys are particle ids and mapped types are vectors of particle ids.
 * @param[in] ghost_contact_pair_candidates Ankerl unordered dense map. Stores
 * potential pairs of local-ghost particle in contact. Keys are particle ids and
 * mapped types are vectors of particle ids.
 * @param[in] sparse_contacts_object The object that contains the
 * information about the mobility status of cells
 */
template <int dim, typename PropertiesIndex>
void
find_particle_particle_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_neighbor_list,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

/**
 * @brief Finds vectors of pairs (particle_particle_candidates) which contains the
 * candidate particle-particle collision pairs. These collision pairs will be
 * used in the fine search to investigate if they are in contact or not.
 *
 * @param[in] particle_handler The particle handler of particles in the broad
 * search.
 * @param[in] cells_local_periodic_neighbor_list A vector (with size equal to
 * the number of local periodic cells at boundary 0) of vectors. Each sub-vector
 * have a size equal to the number of adjacent local cells of the main cell plus
 * one. The first element of each sub-vector shows the main cell itself.
 * @param[in] cells_ghost_periodic_neighbor_list A vector (with size equal to
 * the number of local periodic cells at boundary 0) of vectors. Each sub-vector
 * have a size equal to the number of adjacent ghost cells of the main cell plus
 * one. The first element of each sub-vector shows the main cell itself.
 * @param[in] cells_ghost_local_periodic_neighbor_list A vector (with size equal
 * to the number of ghost periodic cells at boundary 0) of vectors. Each
 * sub-vector have a size equal to the number of adjacent local cells of the
 * main ghost cell plus one. The first element of each sub-vector shows the main
 * ghost cell itself.
 * @param[out] local_contact_pair_periodic_candidates Ankerl unordered dense
 * map. Stores potential pairs of local-local periodic particles in contact.
 * Keys are local particle ids at boundary 0 and mapped types are vectors of
 * local particle ids at boundary 1.
 * @param[out] ghost_contact_pair_periodic_candidates Ankerl unordered dense
 * map. Stores potential pairs of local-ghost periodic particles in contact.
 * Keys are local particle ids at boundary 0 and mapped types are vectors of
 * ghost particle ids at boundary 1.
 * @param[out] ghost_local_contact_pair_periodic_candidates Ankerl unordered
 * dense map. Stores potential pairs of local-ghost periodic particles in
 * contact. Keys are ghost at boundary 0 particle ids and mapped types are
 * vectors of local particle ids at boundary 1.
 */
template <int dim>
void
find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates);

/**
 * @brief Finds a vector of pairs (particle_particle_candidates) which contains the
 * candidate particle-particle collision pairs. These collision pairs will be
 * used in the fine search to investigate if they are in contact or not.
 * This version of the function is used when adaptive sparse contacts is
 * enabled.
 *
 * @param[in] particle_handler The particle handler of particles in the broad
 * search.
 * @param[in] cells_local_periodic_neighbor_list A vector (with size equal to
 * the number of local periodic cells at boundary 0) of vectors. Each sub-vector
 * have a size equal to the number of adjacent local cells of the main cell plus
 * one. The first element of each sub-vector shows the main cell itself.
 * @param[in] cells_ghost_periodic_neighbor_list A vector (with size equal to
 * the number of local periodic cells at boundary 0) of vectors. Each sub-vector
 * have a size equal to the number of adjacent ghost cells of the main cell plus
 * one. The first element of each sub-vector shows the main cell itself.
 * @param[in] cells_ghost_local_periodic_neighbor_list A vector (with size equal
 * to the number of ghost periodic cells at boundary 0) of vectors. Each
 * sub-vector have a size equal to the number of adjacent local cells of the
 * main ghost cell plus one. The first element of each sub-vector shows the main
 * ghost cell itself.
 * @param[out] local_contact_pair_periodic_candidates Ankerl unordered dense
 * map. Stores potential pairs of local-local periodic particles in contact.
 * Keys are local particle ids at boundary 0 and mapped types are vectors of
 * local particle ids at boundary 1.
 * @param[out] ghost_contact_pair_periodic_candidates Ankerl unordered dense
 * map. Stores potential pairs of local-ghost periodic particles in contact.
 * Keys are local particle ids at boundary 0 and mapped types are vectors of
 * ghost particle ids at boundary 1.
 * @param[out] ghost_local_contact_pair_periodic_candidates Ankerl unordered
 * dense map. Stores potential pairs of local-ghost periodic particles in
 * contact. Keys are ghost at boundary 0 particle ids and mapped types are
 * vectors of local particle ids at boundary 1.
 * @param sparse_contacts_object The object that contains the
 * information about the mobility status of cells
 */
template <int dim, typename PropertiesIndex>
void
find_particle_particle_periodic_contact_pairs(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_local_periodic_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_periodic_neighbor_list,
  const typename DEM::dem_data_structures<dim>::cells_neighbor_list
    &cells_ghost_local_periodic_neighbor_list,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &local_contact_pair_periodic_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_contact_pair_periodic_candidates,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &ghost_local_contact_pair_periodic_candidates,
  const AdaptiveSparseContacts<dim, PropertiesIndex> &sparse_contacts_object);

/**
 * @brief Stores the candidate particle-particle collision pairs with a given
 * particle iterator. particle_begin iterator is useful to skip storage of the
 * first particle in main cell (particle_begin will be the iterator after the
 * particles_to_evaluate.begin() in that case). When particle_begin is
 * particles_to_evaluate.begin(), it stores all the particle id in
 * contact_pair_candidates.
 *
 * @param main_particle_id The id of the main particle to store the candidate.
 * @param particle_begin The particle iterator to start storing particle ids.
 * @param particles_to_evaluate The particle range to evaluate with the
 * particle iterator prior storing ids.
 * @param contact_pair_candidates A map which will contain all the particle
 * pairs candidate.
 */
template <int dim>
inline void
store_candidates(
  const types::particle_index &main_particle_id,
  const typename dealii::Particles::ParticleHandler<
    dim>::particle_iterator_range::iterator &particle_begin,
  const typename dealii::Particles::ParticleHandler<
    dim>::particle_iterator_range &particles_to_evaluate,
  typename DEM::dem_data_structures<dim>::particle_particle_candidates
    &contact_pair_candidates);

#endif
