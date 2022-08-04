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

  void
  execute_particle_displacement(
    const Particles::ParticleHandler<dim> &particle_handler);



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
  check_and_move_particles(
    boundary_cells_info_struct<dim> &cell_1,
    boundary_cells_info_struct<dim> &cell_2,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
      &particles_in_cell);

  std::map<types::global_cell_index, types::global_cell_index>
    global_periodic_cell_pair;

  unsigned int direction;
};

#endif /* particle_wall_periodic_displacement_h */