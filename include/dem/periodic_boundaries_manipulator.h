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
 * This class is used to manipulate the particle location dealing with
 * periodic boundaries.
 *
 * @note Currently not fully implemented to work with particle collision, only
 *       allow the periodic cells mapping and particle displacement.
 *
 * @author Audrey Collard-Daigneault, Polytechnique Montreal 2022-
 */

template <int dim>
class PeriodicBoundariesManipulator
{
public:
  PeriodicBoundariesManipulator<dim>();

  /**
   * @brief Set the periodic boundaries parameters. Parameters are implemented
   * to allow use of more than one PBC, but the feature is not implemented yet.
   *
   * @param outlet_boundaries Vector of periodic boundaries identified as outlet
   * @param periodic_boundaries Vector of periodic boundaries identified as
   *                            periodic
   * @param periodic_directions Vector of directions, perpendicular axis of PB
   */
  void
  set_periodic_boundaries_information(
    std::vector<unsigned int> outlet_boundaries,
    std::vector<unsigned int> periodic_boundaries,
    std::vector<unsigned int> periodic_directions)
  {
    outlet_boundary_ids   = outlet_boundaries;
    periodic_boundary_ids = periodic_boundaries;
    directions            = periodic_directions;
  }

  /**
   * @brief Set the periodic boundaries parameters. Parameters are implemented
   * to allow use of more than one PBC, but the feature is not implemented yet
   *
   * @param triangulation Triangulation of mesh
   */
  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * @brief Move particles passing through periodic boundaries (any side)
   * particle_handle doesn't allow automated particle displacement since it is
   * not linked to triangulation and its periodic mapping.
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   */
  void
  execute_particles_displacement(
    const Particles::ParticleHandler<dim> &particle_handler);

private:
  /**
   * @brief Get boundary information related to the face at outlet or periodic
   * boundary and store in boundary_cells_info_struct object.
   *
   * @param cell Current cell on boundary
   * @param face_id Face located on boundary
   * @param boundary_information Reference to the object with boundary info
   */
  void
  get_boundary_info(typename Triangulation<dim>::cell_iterator cell,
                    unsigned int                               face_id,
                    boundary_cells_info_struct<dim> &boundary_information);

  /**
   * @brief Check if particle is outside of the cell, if so, modify the
   * location of the particle with the distance between the periodic faces.
   *
   * @param cell_1 Cell where particles may get outside of domain
   * @param cell_2 Periodic cell to the cell_1
   * @param particles_in_cell Iterator to the particle in cell_1
   */
  void
  check_and_move_particles(
    boundary_cells_info_struct<dim> &cell_1,
    boundary_cells_info_struct<dim> &cell_2,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
      &particles_in_cell);

  std::vector<unsigned int> outlet_boundary_ids;
  std::vector<unsigned int> periodic_boundary_ids;
  std::vector<unsigned int> directions;

  // Mapping the cell pair, cell index is outlet
  std::map<
    types::global_cell_index,
    std::pair<boundary_cells_info_struct<dim>, boundary_cells_info_struct<dim>>>
    periodic_boundary_cells_information;
};

#endif /* particle_wall_periodic_displacement_h */