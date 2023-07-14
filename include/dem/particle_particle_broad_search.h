/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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
 * Author: Shahab Golshan, Polytechnique Montreal, 2019
 */
#include <dem/data_containers.h>
#include <dem/disable_contacts.h>
#include <dem/find_boundary_cells_information.h>
#include <dem/particle_particle_contact_info.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

template <int dim>
class DEMContactManager;

using namespace dealii;

#ifndef particle_particle_broad_search_h
#  define particle_particle_broad_search_h

/**
 * This class is used for broad particle-particle contact search. Broad search
 * is used to obtain all the particle-particle pairs in adjacent cells.
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleParticleBroadSearch
{
public:
  ParticleParticleBroadSearch<dim>();

  /**
   * @brief Finds a vector of pairs (particle_particle_candidates) which shows the
   * candidate particle-particle collision pairs. These collision pairs will be
   * used in the fine search to investigate if they are in contact or not.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair candidates with other
   * containers with neighbors lists
   */

  void
  find_particle_particle_contact_pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    DEMContactManager<dim> &                 container_manager);

  /**
   * @brief Finds a vector of pairs (particle_particle_candidates) which shows the
   * candidate particle-particle collision pairs. These collision pairs will be
   * used in the fine search to investigate if they are in contact or not.
   * This version of the function is used when disabling particle contacts
   * regards mobility is enable.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair candidates with other
   * containers with neighbors lists
   * @param disable_contacts_object The object that contains the
   * information about the mobility status of cells
   */

  void
  find_particle_particle_contact_pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    DEMContactManager<dim> &                 container_manager,
    const DisableContacts<dim> &             disable_contacts_object);

  /**
   * @brief Finds a vector of pairs (particle_particle_candidates) which contains the
   * candidate particle-particle collision pairs. These collision pairs will be
   * used in the fine search to investigate if they are in contact or not.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair periodic candidates with other
   * containers with periodic neighbors lists
   */

  void
  find_particle_particle_periodic_contact_pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    DEMContactManager<dim> &                 container_manager);

  /**
   * @brief Finds a vector of pairs (particle_particle_candidates) which contains the
   * candidate particle-particle collision pairs. These collision pairs will be
   * used in the fine search to investigate if they are in contact or not.
   * This version of the function is used when disabling particle contacts
   * regards mobility is enable.
   *
   * @param particle_handler The particle handler of particles in the broad
   * search
   * @param container_manager The container manager object that contains
   * containers to modify of contact pair periodic candidates with other
   * containers with periodic neighbors lists
   * @param disable_contacts_object The object that contains the
   * information about the mobility status of cells
   */

  void
  find_particle_particle_periodic_contact_pairs(
    dealii::Particles::ParticleHandler<dim> &particle_handler,
    DEMContactManager<dim> &                 container_manager,
    const DisableContacts<dim> &             disable_contacts_object);

private:
  /**
   * @brief Stores the candidate particle-particle collision pairs with a given
   * particle iterator. particle_begin iterator is useful to skip storage of the
   * first particle in main cell (particle_begin will be the iterator after the
   * particles_to_evaluate.begin() in that case). When particle_begin is
   * particles_to_evaluate.begin(), it stores all the particle id in
   * contact_pair_candidates_container.
   *
   * @param particle_begin The particle iterator to start storing particle id
   * @param particles_to_evaluate The particle range to evaluate with the
   * particle iterator prior storing ids
   * @param contact_pair_candidates_container A vector which will contain all
   * the particle pairs which are collision candidates
   */
  inline void
  store_candidates(
    const typename Particles::ParticleHandler<
      dim>::particle_iterator_range::iterator &particle_begin,
    const typename Particles::ParticleHandler<dim>::particle_iterator_range
      &                                 particles_to_evaluate,
    std::vector<types::particle_index> &contact_pair_candidates_container)
  {
    // Create a arbitrary temporary empty container
    if (contact_pair_candidates_container.empty())
      {
        contact_pair_candidates_container.reserve(40);
      }

    // Store particle ids from the selected particle iterator
    for (auto particle_iterator = particle_begin;
         particle_iterator != particles_to_evaluate.end();
         ++particle_iterator)
      {
        contact_pair_candidates_container.emplace_back(
          particle_iterator->get_id());
      }
  }
};

#endif /* particle_particle_broad_search_h */
