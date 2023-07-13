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
#include <dem/data_containers.h>
#include <dem/dem_contact_manager.h>
#include <dem/dem_solver_parameters.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/particle_iterator.h>

#include <iostream>
#include <vector>

using namespace dealii;

#ifndef periodic_boundaries_manipulator_h
#  define periodic_boundaries_manipulator_h

/**
 * This class corresponds to a manipulator of the particles crossing periodic
 * cells in DEM. It maps cell information of the pairs of periodic cells
 * and manipulate the particles location when they cross periodic boundaries,
 * i.e. when a particle cross a periodic boundary (0 or 1) its location is
 * modified with the offset between faces of periodic cells on periodic
 * boundaries.
 *
 * @note Currently the code only supports one pair of periodic boundaries and
 * those boundaries must be parallel and aligned with an axis of the domain
 * (e.g., x axis).
 */

template <int dim>
class PeriodicBoundariesManipulator
{
public:
  PeriodicBoundariesManipulator<dim>();

  /**
   * @brief Sets the periodic boundaries parameters
   *
   * @param periodic_boundary_id_0 Id of the first periodic boundary
   * @param periodic_boundary_id_1 Id of the second periodic boundary
   * @param periodic_direction Perpendicular axis of PB
   */
  void
  set_periodic_boundaries_information(
    const types::boundary_id periodic_boundary_id_0,
    const types::boundary_id periodic_boundary_id_1,
    const unsigned int       periodic_direction)
  {
    periodic_boundary_0 = periodic_boundary_id_0;
    periodic_boundary_1 = periodic_boundary_id_1;
    direction           = periodic_direction;
  }

  /**
   * @brief Execute the mapping of the cells on periodic boundaries and store
   * information in periodic_boundaries_cells_information
   *
   * @param triangulation Triangulation of mesh
   * @param periodic_boundaries_cells_information Map of information of the
   * pair of cells on periodic boundaries
   */
  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information);

  /**
   * @brief Moves particles crossing periodic boundaries (any side)
   * particle_handle doesn't allow automated particle displacement since it is
   * not linked to triangulation and its periodic mapping
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   * @param periodic_boundaries_cells_information Map of information of the
   * pair of cells on periodic boundaries
   */
  void
  execute_particles_displacement(
    const Particles::ParticleHandler<dim> &particle_handler,
    typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information);

  /**
   * @brief Return the periodic offset distance, it is calculated from the first pair of
   * cells on periodic boundaries, all pair of cells are assumed to have the
   * same offset
   */
  inline Tensor<1, dim>
  get_periodic_offset_distance()
  {
    return constant_periodic_offset;
  }

private:
  /**
   * @brief Gets boundaries information related to the face at periodic
   * boundaries 0 and 1 and stores in periodic_boundaries_cells_info_struct
   * object
   *
   * @param cell Current cell on boundary
   * @param face_id Face located on boundary
   * @param boundaries_information Reference to the object with periodic
   * boundaries information
   */
  void
  get_periodic_boundaries_info(
    typename Triangulation<dim>::cell_iterator  cell,
    const unsigned int                          face_id,
    periodic_boundaries_cells_info_struct<dim> &boundaries_information);

  /**
   * @brief Checks if particle is outside of the cell, if so, modifies the
   * location of the particle with the distance (offset) between the periodic
   * faces
   *
   * @param boundaries_cells_content Reference to the object with periodic
   * boundaries information
   * @param particles_in_cell Iterator to the particles in cell
   * @param particles_in_pb0_cell If the particles are linked to a cell on the
   * periodic boundary 0 (true) or the periodic boundary 1 (false)
   */
  void
  check_and_move_particles(
    const periodic_boundaries_cells_info_struct<dim> &boundaries_cells_content,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
      &        particles_in_cell,
    const bool particles_in_pb0_cell);

  types::boundary_id periodic_boundary_0;
  types::boundary_id periodic_boundary_1;
  unsigned int       direction;
  Tensor<1, dim>     constant_periodic_offset;
};



#endif /* periodic_boundaries_manipulator_h */
