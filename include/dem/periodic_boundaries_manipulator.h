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
 * This class is used to manipulate the particles location when they cross a
 * periodic boundary
 *
 * @note Particle collisions across periodic boundaries are currently not
 * implemented.

 */

template <int dim>
class PeriodicBoundariesManipulator
{
public:
  PeriodicBoundariesManipulator<dim>();

  /**
   * @brief Sets the periodic boundaries parameters
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
   * @brief Sets the periodic boundaries parameters
   *
   * @param triangulation Triangulation of mesh
   */
  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation);

  /**
   * @brief Moves particles crossing periodic boundaries (any side)
   * particle_handle doesn't allow automated particle displacement since it is
   * not linked to triangulation and its periodic mapping
   *
   * @param particle_handler Particle handler of particles located in boundary
   * cells
   */
  void
  execute_particles_displacement(
    const Particles::ParticleHandler<dim> &particle_handler);

private:
  /**
   * @brief Gets boundaries information related to the face at outlet and
   * periodic boundaries and stores in  periodic_boundaries_cells_info_struct
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
   * location of the particle with the distance between the periodic faces.
   *
   * @param boundaries_cells_content Reference to the object with periodic
   * boundaries information
   * @param particles_in_cell Iterator to the particles in cell
   * @param particles_in_outlet_cell If the particles are linked to the outlet
   * cell or the periodic cell
   */
  void
  check_and_move_particles(
    const periodic_boundaries_cells_info_struct<dim> &boundaries_cells_content,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
              &particles_in_cell,
    const bool particles_in_outlet_cell);

  std::vector<unsigned int> outlet_boundary_ids;
  std::vector<unsigned int> periodic_boundary_ids;
  std::vector<unsigned int> directions;

  // Mapping of periodic boundaries information, cell index at outlet is the
  // map key
  typename dem_data_containers::dem_data_structures<
    dim>::periodic_boundaries_cells_info periodic_boundaries_cells_information;
};

#endif /* particle_wall_periodic_displacement_h */
