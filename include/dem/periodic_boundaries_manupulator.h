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
 * Author: Audrey Collard-Daigneault, Polytechnique Montreal, 2022
 */

#include <dem/boundary_cells_info_struct.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef particle_wall_periodic_displacement_h
#  define particle_wall_periodic_displacement_h

/**
 * This class is used for ...
 *
 * @note
 *
 * @author Audrey Collard-Daigneault, Polytechnique Montreal 2022-
 */

template <int dim>
class PeriodicBoundariesManipulator
{
public:
  PeriodicBoundariesManipulator<dim>();

  void
  set_periodic_boundaries_direction(unsigned int direction)
  {
    this->direction = direction;
  }

  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation);

  /* void
   find_cell_neighbors_and_periodic(
     const parallel::distributed::Triangulation<dim> &triangulation,
     std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
       &cells_local_periodic_neighbor_list,
     std::vector<std::vector<typename Triangulation<dim>::active_cell_iterator>>
       &cells_ghost_periodic_neighbor_list);*/

  void
  execute_particle_displacement(
    const Particles::ParticleHandler<dim> &particle_handler);



  /* void
   find_particle_particle_periodic_contact_pairs(
     dealii::Particles::ParticleHandler<dim> &particle_handler,
     const std::map<
       int,
       std::pair<boundary_cells_info_struct<dim>,
 boundary_cells_info_struct<dim>>> &periodic_boundaries_info, const std::vector<
       std::vector<typename Triangulation<dim>::active_cell_iterator>>
       *cells_local_periodic_neighbor_list,
     const std::vector<
       std::vector<typename Triangulation<dim>::active_cell_iterator>>
       *cells_ghost_periodic_neighbor_list,
     std::unordered_map<types::particle_index,
                        std::vector<std::pair<types::particle_index, Point<3>>>>
       &local_periodic_contact_pair_candidates,
     std::unordered_map<types::particle_index,
                        std::vector<std::pair<types::particle_index, Point<3>>>>
       &ghost_periodic_contact_pair_candidates); */

private:
  std::map<
    types::global_cell_index,
    std::pair<boundary_cells_info_struct<dim>, boundary_cells_info_struct<dim>>>
    periodic_boundary_cells_information;

  void
  get_boundary_info(typename Triangulation<dim>::cell_iterator cell,
                    unsigned int                               face_id,
                    boundary_cells_info_struct<dim> &boundary_information);

  void
  move_particles(
    boundary_cells_info_struct<dim> &cell_1,
    boundary_cells_info_struct<dim> &cell_2,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
      &particles_in_cell);

  std::map<types::global_cell_index, types::global_cell_index>
    periodic_cell_pair;

  double       translation_value;
  unsigned int direction;
};

#endif /* particle_wall_periodic_displacement_h */