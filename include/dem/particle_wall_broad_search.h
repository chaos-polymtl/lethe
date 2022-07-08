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

#include <dem/boundary_cells_info_struct.h>
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
 * This class is used for broad particle-wall contact search. Broad search
 * is used to obtain all the particles located at boundary cells
 *
 * @note
 *
 * @author Shahab Golshan, Polytechnique Montreal 2019-
 */

template <int dim>
class ParticleWallBroadSearch
{
public:
  ParticleWallBroadSearch<dim>();

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
      &                                    boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    std::unordered_map<
      types::particle_index,
      std::unordered_map<unsigned int,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>,
                                    types::boundary_id,
                                    unsigned int>>>
      &particle_wall_contact_candidates);

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
   * @param pfw_contact_candidates Output of particle-floating wall broad search
   * which contains all the particle-floating wall collision candidates
   */

  void
  find_particle_floating_wall_contact_pairs(
    const std::unordered_map<
      types::particle_index,
      std::set<typename Triangulation<dim>::active_cell_iterator>>
      &                                    boundary_cells_for_floating_walls,
    const Particles::ParticleHandler<dim> &particle_handler,
    const Parameters::Lagrangian::FloatingWalls<dim> &floating_wall_properties,
    const double &                                    simulation_time,
    std::unordered_map<
      types::particle_index,
      std::unordered_map<unsigned int, Particles::ParticleIterator<dim>>>
      &pfw_contact_candidates);
};


/**
 * UPDATE THIS *****************************
 *
 * @note
 *
 */

template <int dim>
class ParticleMovingMeshBroadSearch
{
public:
  ParticleMovingMeshBroadSearch<dim>();

  /**
   * UPDATE ****************************
   *
   * @param
   */

  void
  find_particle_moving_mesh_contact_pairs(
    const std::unordered_map<
      typename Triangulation<dim>::active_cell_iterator,
      std::unordered_map<int,
                         typename Triangulation<dim>::active_cell_iterator>>
      &                                    moving_mesh_information,
    const Particles::ParticleHandler<dim> &particle_handler,
    std::unordered_map<
      unsigned int,
      std::unordered_map<types::particle_index,
                         std::tuple<Particles::ParticleIterator<dim>,
                                    Tensor<1, dim>,
                                    Point<dim>>>>
      &particle_moving_mesh_contact_candidates);
};

#endif /* particle_wall_broad_search_h */
