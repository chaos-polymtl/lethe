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

#include <core/data_containers.h>

#include <dem/dem_properties.h>
#include <dem/dem_solver_parameters.h>
#include <dem/particle_wall_contact_info_struct.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_wall_fine_search_h
#  define particle_wall_fine_search_h

/**
 * This class is used for fine particle-wall contact search. Fine search
 * is used to find all the particles which are physically in contact with
 * system boundaries and obtain all the required information for calculation
 * of the corresponding particle-wall contact force
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleWallFineSearch
{
public:
  ParticleWallFineSearch<dim>();

  /**
   * Iterates over the contact candidates from particle-wall broad search
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

  void
  particle_wall_fine_search(
    const typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_candidates &particle_wall_contact_pair_candidates,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_in_contact &particle_wall_pairs_in_contact);

  /**
   * Iterates over the contact candidates from particle-floating wall broad
   * search (pfw_contact_candidates) to add new contact pairs to the
   * pfw_pairs_in_contact container
   *
   * @param pfw_contact_pair_candidates The output of particle-floating wall
   * broad search which shows contact pair candidates
   * @param floating_wall_properties Properties of floating walls defined in the
   * parameter handler
   * @param simulation_time Simulation time
   * @param pfw_pairs_in_contact An unordered_map of maps which stores
   * all the particle-floating wall pairs which are in contact, and
   * the contact information in a struct. Note that the size of this unordered
   * map is equal to the number of particles
   */
  void
  particle_floating_wall_fine_search(
    const typename dem_data_containers::dem_data_structures<
      dim>::particle_floating_wall_candidates &       pfw_contact_candidates,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double &                                    simulation_time,
    typename dem_data_containers::dem_data_structures<
      dim>::particle_wall_in_contact &pfw_pairs_in_contact);


  /**
   * Iterates over the contact candidates from particle-floating mesh broad
   * search (particle_floating_mesh_contact_candidates) to add new contact pairs
   * to the particle_floating_mesh_in_contact container
   *
   * @param particle_floating_mesh_contact_candidates The output of particle-floating mesh
   * broad search which shows contact pair candidates
   * @param particle_floating_mesh_in_contact An map of maps which stores
   * all the particle-floating mesh pairs which are in contact
   */

  void
  particle_floating_mesh_fine_search(
    const typename dem_data_containers::dem_data_structures<
      dim>::particle_floating_mesh_candidates
      &particle_floating_mesh_contact_candidates,
    typename dem_data_containers::dem_data_structures<dim>::
      particle_floating_mesh_in_contact &particle_floating_mesh_in_contact);
};

#endif /* particle_wall_fine_search_h */
