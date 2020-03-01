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
#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>
#include <deal.II/particles/property_pool.h>

#include <dem/boundary_cells_info_struct.h>
#include <dem/dem_solver_parameters.h>
#include <dem/integrator.h>
#include <dem/nonuniform_insertion.h>
#include <dem/pp_broad_search.h>
#include <dem/pp_contact_force.h>
#include <dem/pp_contact_info_struct.h>
#include <dem/pp_fine_search.h>
#include <dem/pw_broad_search.h>
#include <dem/pw_contact_force.h>
#include <dem/pw_contact_info_struct.h>
#include <dem/pw_fine_search.h>
#include <dem/uniform_insertion.h>
#include <dem/visualization.h>
#include <dem/write_vtu.h>

#include <chrono>
#include <math.h>

#ifndef DEMITERATOR_H_
#define DEMITERATOR_H_

/**
 * Iterator on the functions and classes for each DEM step
 */

template <int dim> class DEM_iterator {
public:
  DEM_iterator<dim>();

  /**
   * DEM iterator which handles all the operations at each DEM step
   *
   * @param particle_handler Particle handler to access all the particles in the
   * system
   * @param triangulation Triangulation of the DEM system
   * @param DEM_step Current DEM step
   * @param DEM_time Current DEM time
   * @param cell_neighbor_list A vector which contains all the neighbors of
   * cells present in the triangulation
   * @param pairs_in_contact_info A vector of maps which contains information of
   * particle pairs in contact
   * @param boundary_cells_information A vector which contains the information
   * of boundary cells and faces
   * @param pw_pairs_in_contact A vector of maps which contains information of
   * particle-wall pairs in contact
   * @param dem_parameters DEM parameters
   * @param g Gravitational acceleration
   * @param properties Properties of particles which are used in the
   * visualization
   * @param property_pool Property pool of particles in the system
   * @param pp_contact_force_object Particle-particle contact force object
   * @param pw_contact_force_object Particle-wall contact force object
   * @param integrator_object Integration object
   * @param pp_broad_search_object Particle-particle broad search object
   * @param pp_fine_search_object Particle-particle fine search object
   * @param pw_broad_search_object Particle-wall broad search object
   * @param pw_fine_search_object Particle-wall fine search object
   */
  void engine(
      Particles::ParticleHandler<dim> &particle_handler,
      const Triangulation<dim> &triangulation, int &DEM_step, double &DEM_time,
      std::vector<std::set<typename Triangulation<dim>::active_cell_iterator>>
          &cell_neighbor_list,
      std::vector<std::map<int, pp_contact_info_struct<dim>>>
          &pairs_in_contact_info,
      std::vector<boundary_cells_info_struct<dim>> &boundary_cells_information,
      std::vector<std::map<int, pw_contact_info_struct<dim>>>
          &pw_pairs_in_contact,
      DEMSolverParameters<dim> &dem_parameters, Tensor<1, dim> &g,
      std::vector<std::pair<std::string, int>> properties,
      Particles::PropertyPool &property_pool,
      PPContactForce<dim> *pp_contact_force_object,
      PWContactForce<dim> *pw_contact_force_object,
      Integrator<dim> *integrator_object,
      PPBroadSearch<dim> *pp_broad_search_object,
      PPFineSearch<dim> *pp_fine_search_object,
      PWBroadSearch<dim> *pw_broad_search_object,
      PWFineSearch<dim> *pw_fine_search_object);

private:
  /**
   * Reinitializes exerted forces and momentums on particles
   *
   * @param particle_handler Particle handler to access all the particles in the
   * system
   */
  void reinitialize_force(Particles::ParticleHandler<dim> &particle_handler);
};

#endif /* DEMITERATOR_H_ */
