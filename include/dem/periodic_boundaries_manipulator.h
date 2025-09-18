// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundaries_manipulator_h
#define lethe_periodic_boundaries_manipulator_h

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>


using namespace dealii;

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
  PeriodicBoundariesManipulator();

  /**
   * @brief Sets the periodic boundaries parameters.
   *
   * @param[in] periodic_boundary_id_0 ID of the first periodic boundary.
   * @param[in] periodic_direction Perpendicular axis of PB.
   */
  void
  set_periodic_boundaries_information(
    const types::boundary_id periodic_boundary_id_0,
    const unsigned int       periodic_direction)
  {
    // If function is reached, there are periodic boundaries in the simulation
    periodic_boundaries_enabled = true;

    // Communicate to the action manager that there are periodic boundaries
    DEMActionManager::get_action_manager()->set_periodic_boundaries_enabled();

    periodic_boundary_0 = periodic_boundary_id_0;
    direction           = periodic_direction;
  }

  /**
   * @brief Execute the mapping of the cells on periodic boundaries and store
   * information in periodic_boundaries_cells_information.
   *
   * @param[in] triangulation Triangulation.
   * @param[out] periodic_boundaries_cells_information Map of information of the
   * pair of cells on periodic boundaries.
   */
  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information);

  /**
   * @brief Moves particles crossing periodic boundaries (any side)
   * particle_handler doesn't allow automated particle displacement since it is
   * not linked to triangulation and its periodic mapping.
   *
   * @param[in,out] particle_handler Particle handler of particles located in
   * boundary cells.
   * @param[in] periodic_boundaries_cells_information Map of information of the
   * pair of cells on periodic boundaries.
   *
   * @return Flag if at least one particle has been moved to the other periodic
   * cell.
   */
  bool
  execute_particles_displacement(
    const Particles::ParticleHandler<dim> &particle_handler,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information);

  /**
   * @brief Return the periodic offset distance, it is calculated from the first
   * pair of cells on periodic boundaries, all pair of cells are assumed to
   * have the same offset.
   *
   * @return Offset distance between periodic boundaries.
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
   * object.
   *
   * @param[in] cell Current cell on boundary.
   * @param[in] face_id Face located on boundary.
   * @param[out] boundaries_information Reference to the object with periodic
   * boundary information.
   */
  void
  get_periodic_boundaries_info(
    typename Triangulation<dim>::cell_iterator  cell,
    const unsigned int                          face_id,
    periodic_boundaries_cells_info_struct<dim> &boundaries_information);

  /**
   * @brief Checks if particle is outside the cell, if so, modifies the
   * location of the particle with the distance (offset) between the periodic
   * faces.
   *
   * @param[in] boundaries_cells_content Reference to the object with periodic
   * boundary information.
   * @param[in] particles_in_pb0_cell If the particles are linked to a cell on
   * the periodic boundary 0 (true) or the periodic boundary 1 (false).
   * @param[in,out] particles_in_cell Iterator to the particles in cell.
   * @param[out] particle_has_been_moved Flag the particle has been moved to
   * the other periodic cell.
   */
  void
  check_and_move_particles(
    const periodic_boundaries_cells_info_struct<dim> &boundaries_cells_content,
    const bool                                       &particles_in_pb0_cell,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
         &particles_in_cell,
    bool &particle_has_been_moved);

  /**
   * @brief Flag for periodic boundary conditions in simulation. Useful to
   * exit function when there are no periodic boundaries.
   */
  bool periodic_boundaries_enabled;

  /**
   * @brief ID of the first periodic boundary. No needs to store the second one
   * since there are linked on the triangulation, and accessible through
   * functions on cells on the boundary condition 0.
   */
  types::boundary_id periodic_boundary_0;

  /**
   * @brief Direction of the periodic boundary, it is the perpendicular axis of
   * the periodic boundaries.
   */
  unsigned int direction;

  /**
   * @brief Offset distance between periodic boundaries, it is calculated from
   * the first pair of cells on periodic boundaries, all pair of cells are
   * assumed to have the same offset.
   */
  Tensor<1, dim> constant_periodic_offset;
};

#endif
