// SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundaries_manipulator_h
#define lethe_periodic_boundaries_manipulator_h

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>
#include <unordered_map>

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
   * @param[in] periodic_boundary_ids_0 Map of IDs of the first boundary of each PB pair.
   * @param[in] periodic_directions Map of perpendicular axes of each PB pair.
   */
  void
  set_periodic_boundaries_information(
<<<<<<< HEAD
    const std::unordered_map<unsigned int, types::boundary_id> &periodic_boundary_ids_0,
    const std::unordered_map<unsigned int, unsigned int>       &periodic_directions)
=======
    const std::unordered_map<unsigned int, types::boundary_id> periodic_boundary_id_0,
    const std::unordered_map<unsigned int, unsigned int>       periodic_direction)
>>>>>>> bd6a750b8 (merge avec le reste de la branche)
  {
    // If function is reached and vectors are not empty
    if (!periodic_boundary_ids_0.empty())
      {
        periodic_boundaries_enabled = true;

        // Communicate to the action manager that there are periodic boundaries
        DEMActionManager::get_action_manager()->set_periodic_boundaries_enabled();

        periodic_boundaries_ids = periodic_boundary_ids_0;
        directions              = periodic_directions;

        // Initialize offset map
        periodic_offsets.clear();
      }
  }

  /**
   * @brief Execute the mapping of the cells on periodic boundaries and store
   * information in periodic_boundaries_cells_information.
   *
   * @param[in] triangulation Triangulation.
   * @param[out] periodic_boundaries_cells_information Map (multimap) of information of the
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
   * @brief Return the periodic offset distance for a specific boundary ID.
   * It is calculated from the first pair of cells on periodic boundaries, 
   * all pair of cells are assumed to have the same offset.
   *
   * @param[in] boundary_id ID of the boundary to query
   * @return Offset distance between periodic boundaries.
   */
  inline Tensor<1, dim>
  get_periodic_offset_distance(const types::boundary_id boundary_id)
  {
    if (periodic_offsets.find(boundary_id) != periodic_offsets.end())
      return periodic_offsets[boundary_id];
    else
      return Tensor<1, dim>();
  }

private:
  /**
   * @brief Gets boundaries information related to the face at periodic
   * boundaries and stores in periodic_boundaries_cells_info_struct
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
   * location of the particle using the specific offset for boundary.
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
   * @brief IDs of the first periodic boundary of each pair.
   */
<<<<<<< HEAD
  std::unordered_map<unsigned int, types::boundary_id> periodic_boundaries_ids;
=======
  std::unordered_map<unsigned int, types::boundary_id> periodic_boundary_0;
>>>>>>> bd6a750b8 (merge avec le reste de la branche)

  /**
   * @brief Direction of the periodic boundaries, it is the perpendicular axis of
   * the periodic boundaries.
   */
<<<<<<< HEAD
  std::unordered_map<unsigned int, unsigned int> directions;
=======
  std::unordered_map<unsigned int, unsigned int> direction;
>>>>>>> bd6a750b8 (merge avec le reste de la branche)

  /**
   * @brief Map storing offset distance between periodic boundaries, keyed by the
   * boundary ID (pb0). It is calculated from the first pair of cells on periodic
   * boundaries, all pair of cells are assumed to have the same offset.
   */
  std::unordered_map<types::boundary_id, Tensor<1, dim>> periodic_offsets;
};

#endif
