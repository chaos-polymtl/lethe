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
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> & triangulation,
    const std::vector<unsigned int> &                 outlet_boundaries);

    void
  execute_particle_displacement(
    const std::map<int,
                   std::pair<boundary_cells_info_struct<dim>,
                             boundary_cells_info_struct<dim>>>
      &                                    periodic_boundary_cells_information,
    const Particles::ParticleHandler<dim> &particle_handler);

private:
  std::map<
    int,
    std::pair<boundary_cells_info_struct<dim>, boundary_cells_info_struct<dim>>>
    periodic_boundary_cells_information;
};

#endif /* particle_wall_periodic_displacement_h */