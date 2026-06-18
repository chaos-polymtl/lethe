// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundaries_manipulator_h
#define lethe_periodic_boundaries_manipulator_h

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#include <map>
#include <unordered_map>
#include <vector>

using namespace dealii;

/**
 * This class corresponds to a manipulator of the particles crossing periodic
 * cells in DEM. It maps cell information of the pairs of periodic cells
 * and manipulate the particles' location when they cross periodic boundaries,
 * i.e. when a particle cross a periodic boundary (0 or 1) its location is
 * modified with the offset between faces of periodic cells on periodic
 * boundaries.
 *
 * @note Currently the code only supports periodic boundaries that are
 * parallel and aligned with an axis of the domain (e.g., x axis).
 */

template <int dim>
class PeriodicBoundariesManipulator
{
public:
  PeriodicBoundariesManipulator();

  /**
   * @brief Sets the periodic boundaries parameters.
   *
   * @param[in] periodic_directions Map of perpendicular axes of each PB pair,
   * keyed by the principal periodic boundary id (periodic id 0). The keys are
   * the principal periodic boundary ids.
   */
  void
  set_periodic_boundaries_information(
    const std::map<types::boundary_id, unsigned int> &periodic_directions)
  {
    // If this function is reached and the map is not empty
    if (periodic_directions.empty())
      return;

    periodic_boundaries_enabled = true;

    // Communicate to the action manager that there are periodic boundaries
    DEMActionManager::get_action_manager()->set_periodic_boundaries_enabled();

    this->directions = periodic_directions;

    // Initialize offset map
    this->periodic_offsets.clear();
  }

  /**
   * @brief Execute the mapping of the cells on periodic boundaries and store
   * information in periodic_boundaries_cells_information.
   *
   * @param[in] triangulation Triangulation.
   * @param[out] periodic_boundaries_cells_information Map (multimap) of
   * information of the pair of cells on periodic boundaries.
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
   * @brief Return the periodic offset distance for a specific pb0 boundary ID.
   * All pairs of cells on that pair of periodic boundaries are assumed to have
   * the same offset. If no periodic offset has been identified for the
   * boundary id, a zero tensor is inserted as a default value.
   *
   * @param[in] boundary_id ID of the boundary to query.
   *
   * @return Offset distance between periodic boundaries.
   */
  inline const Tensor<1, dim> &
  get_periodic_offset_distance(const types::boundary_id boundary_id)
  {
    // If no offset was identified during looping over the cell, a zero
    // tensor is returned as a default offset.
    if (!periodic_offsets.contains(boundary_id))
      {
        Tensor<1, dim> zero_tensor;
        periodic_offsets.insert(std::make_pair(boundary_id, zero_tensor));
      }
    return periodic_offsets.at(boundary_id);
  }

  /**
   * @brief Return the directions of the periodic boundary pairs, keyed by the
   * principal periodic boundary id (periodic id 0). The keys are the principal
   * periodic boundary ids.
   *
   * @return Map of periodic directions keyed by principal periodic boundary id.
   */
  inline const std::map<types::boundary_id, unsigned int> &
  get_periodic_directions() const
  {
    return directions;
  }

  /**
   * @brief Return the combined periodic offsets
   *
   * @return Combined periodic offsets.
   */
  inline const std::vector<Tensor<1, dim>> &
  get_combined_periodic_offsets() const
  {
    return combined_periodic_offsets;
  }

private:
  /**
   * @brief Get boundary information related to a face at a periodic
   * boundary and store it in a periodic_boundaries_cells_info_struct
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
   * @brief Check if particle is outside the cell, if so, modify the
   * location of the particle using the relevant periodic offset for the
   * boundary.
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
   * @brief Compute the combined periodic offsets and store them in
   * combined_periodic_offsets.
   */
  void
  compute_combined_periodic_offsets();

  /**
   * @brief Flag for periodic boundary conditions in simulation. Useful to
   * exit function when there are no periodic boundaries.
   */
  bool periodic_boundaries_enabled;

  /**
   * @brief Direction of the periodic boundaries, it is the perpendicular axis
   * of the periodic boundaries. Keys of this map are the principal periodic
   * boundary ids (periodic id 0).
   */
  std::map<types::boundary_id, unsigned int> directions;

  /**
   * @brief Map storing offset distance between periodic boundaries, keyed by
   * the boundary ID (pb0). It is calculated from the first pair of cells on
   * periodic boundaries, so all pairs of cells on a given peridodic boundary
   * are assumed to have the same offset.
   */
  std::unordered_map<types::boundary_id, Tensor<1, dim>> periodic_offsets;

  /**
   * @brief Storage for all 9 (2D) or 27 (3D) precomputed periodic translation
   * vectors. Calculated from periodic_offsets.
   */
  std::vector<Tensor<1, dim>> combined_periodic_offsets;
};

#endif
