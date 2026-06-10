// SPDX-FileCopyrightText: Copyright (c) 2022-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_periodic_boundaries_manipulator_h
#define lethe_periodic_boundaries_manipulator_h

#include <dem/boundary_cells_info_struct.h>
#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/particles/particle_handler.h>

#include <unordered_map>
#include <vector>

using namespace dealii;

/**
 * This class corresponds to a manipulator of the particles crossing periodic
 * cells in DEM. It maps cell information of the pairs of periodic cells
 * and manipulates the particles' location when they cross periodic boundaries,
 * i.e. when a particle crosses a periodic boundary (0 or 1) its location is
 * modified with the offset between faces of periodic cells on periodic
 * boundaries.
 *
 * @note PBC must be  parallel and aligned with an axis of the domain
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
   * @param[in] periodic_boundary_ids_0 Map of IDs of the primary boundary of
   * each PB pair.
   * @param periodic_boundary_ids_1 Map of IDs of the secondary boundary of
   * each PB pair.
   * @param[in] periodic_directions Map of perpendicular axes of each PB pair.
   * @param[in] prm_boundary_index Index of the PB conditions in the .prm file.
   */
  void
  set_periodic_boundaries_information(
    const std::unordered_map<unsigned int, types::boundary_id>
      &periodic_boundary_ids_0,
    const std::unordered_map<unsigned int, types::boundary_id>
      &periodic_boundary_ids_1,
    const std::unordered_map<unsigned int, unsigned int> &periodic_directions,
    const std::vector<unsigned int>                      &prm_boundary_index)
  {
    // If this function is reached and vectors are not empty
    if (periodic_boundary_ids_0.empty())
      return;

    periodic_boundaries_enabled = true;

    // Communicate to the action manager that there are periodic boundaries
    DEMActionManager::get_action_manager()->set_periodic_boundaries_enabled();

    primary_periodic_boundaries_ids   = periodic_boundary_ids_0;
    secondary_periodic_boundaries_ids = periodic_boundary_ids_1;
    directions                        = periodic_directions;
    prm_periodic_bc_index             = prm_boundary_index;

    number_of_declared_periodic_boundaries = prm_periodic_bc_index.size();
    boundary_combination_set_to_index = build_periodic_boundary_combinations();
  }

  using BoundarySet = std::set<types::boundary_id>;
  /**
   * @brief Build all valid periodic boundary combinations and assign a unique
   * index to each combination.
   *
   * Example, we have two PBC, between boundary 0 to 1 and 6 to 3, respectively:
   *   PBC 0: (0 <-> 1)
   *   PBC 1: (6 <-> 3)
   *
   * Valid combinations generated:
   *
   *   {0}
   *   {1}
   *   {6}
   *   {3}
   *   {0,6}
   *   {0,3}
   *   {1,6}
   *   {1,3}
   *
   * Invalid combinations such as {0,1} are NOT assigned an index because both
   * boundaries belong to the same periodic pair. In other words, a cell will
   * never touch both boundaries 0 and 1. (Except if the triangulation is one
   * cell in thickness, which should not happen)
   *
   * The generated map associates:
   *
   *   BoundarySet -> unique integer index
   */
  std::map<BoundarySet, std::uint8_t>
  build_periodic_boundary_combinations()
  {
    // Store all the boundary periodic boundary condition IDs like they are
    // declared in the prm file.
    std::vector<unsigned int> prm_pbc_ids(prm_periodic_bc_index.begin(),
                                          prm_periodic_bc_index.end());

    // Sort IDs to guarantee deterministic ordering of generated indices.
    // NOTE FROM GABO: I am pretty sure this is not necessary.
    std::ranges::sort(prm_pbc_ids);

    // Initialize the container.
    std::map<BoundarySet, std::uint8_t> combinations_to_id;

    // Temporary set used during recursive construction
    BoundarySet current;

    // Counter used to assign the ID to each combination associated with each
    // valid combination.
    std::uint8_t counter = 0;

    // Lambda function. This function will be called recursively.
    auto build_valid_combination = [&](auto            &&self,
                                       const std::size_t pbc_index) -> void //
    {
      if (pbc_index == prm_pbc_ids.size())
        {
          if (!current.empty())
            {
              // Use emplace to construct the pair in-place
              combinations_to_id.emplace(current, counter++);
            }
          return;
        }

      const auto pbc_id = prm_pbc_ids[pbc_index];

      self(self, pbc_index + 1);

      const auto primary_id = primary_periodic_boundaries_ids.at(pbc_id);
      current.insert(primary_id);
      self(self, pbc_index + 1);
      current.erase(primary_id); // Backtrack

      const auto secondary_id = secondary_periodic_boundaries_ids.at(pbc_id);
      current.insert(secondary_id);
      self(self, pbc_index + 1);
      current.erase(secondary_id); // Backtrack
    };

    // Kick off recursion only if there are elements to process
    build_valid_combination(build_valid_combination, 0);

    return combinations_to_id;
  }

  /**
   * @brief Returns the periodic boundary ID set that is common between two cell
   * iterator.
   *
   * @param main_boundaries All the periodic boundary ID(s) of the main cell.
   * @param neighbor_boundaries All the periodic boundary  ID(s) of the
   * neighboring cell.
   */
  BoundarySet
  find_shared_periodic_boundaries(const BoundarySet &main_boundaries,
                                  const BoundarySet &neighbor_boundaries) const
  {
    BoundarySet shared_pbc_ids;

    for (const auto &[pbc_id, primary_id] : primary_periodic_boundaries_ids)
      {
        const auto secondary_id = secondary_periodic_boundaries_ids.at(pbc_id);

        const bool match = (main_boundaries.contains(primary_id) &&
                            neighbor_boundaries.contains(secondary_id)) ||
                           (main_boundaries.contains(secondary_id) &&
                            neighbor_boundaries.contains(primary_id));

        if (match)
          {
            if (main_boundaries.contains(primary_id))
              shared_pbc_ids.insert(primary_id);
            else
              shared_pbc_ids.insert(secondary_id);
          }
      }

    return shared_pbc_ids;
  }

  /**
   * @brief Retrieves the index of a container based on its identifier.
   *
   * @param share_combination_of_boundaries The boundary id set in common between
   * the main and neighboring cells..
   * @return The index associated with the right set.
   */
  std::uint8_t
  get_container_index(const BoundarySet &share_combination_of_boundaries) const
  {
    return boundary_combination_set_to_index.at(
      share_combination_of_boundaries);
  }

  /**
   * @brief Execute the mapping of the cells on periodic boundaries and store
   * information in periodic_boundaries_cells_information.
   *
   * @param[in] triangulation Triangulation.
   * @param[out] periodic_boundaries_cells_information Map (multimap) of
   * information of the pair of cells on periodic boundaries.
   * @param cell_to_pbc_mesh_id_set
   */
  void
  map_periodic_cells(
    const parallel::distributed::Triangulation<dim> &triangulation,
    typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information,
    typename DEM::dem_data_structures<dim>::cell_touch_boundary_id
      &cell_to_pbc_mesh_id_set);

  /**
   * @brief Moves particles crossing periodic boundaries (any side)
   * particle_handler doesn't allow automated particle displacement since it is
   * not linked to triangulation and its periodic mapping.
   *
   * @param[in,out] particle_handler Particle handler of particles located in
   * boundary cells.
   * @param[in] periodic_boundaries_cells_information Map of information of the
   * pair of cells on periodic boundaries.
   * @param particle_to_total_periodic_displacement
   *
   * @return Flag if at least one particle has been moved to the other periodic
   * cell.
   */
  bool
  execute_particles_displacement(
    const Particles::ParticleHandler<dim> &particle_handler,
    const typename DEM::dem_data_structures<dim>::periodic_boundaries_cells_info
      &periodic_boundaries_cells_information,
    typename DEM::dem_data_structures<dim>::particle_index_tensor_map
      &particle_to_total_periodic_displacement);

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
   * @brief Return the index of the periodic boundary conditions in the .prm.
   *
   * @return Index of the periodic boundary conditions.
   */
  inline const std::vector<unsigned int> &
  get_primary_periodic_bc_index() const
  {
    return prm_periodic_bc_index;
  }

  /**
   * @brief Return the mesh IDs of the principal periodic boundaries
   *
   * @return Mesh IDS of the principal periodic boundaries.
   */
  inline const std::unordered_map<unsigned int, types::boundary_id> &
  get_periodic_boundaries_ids() const
  {
    return primary_periodic_boundaries_ids;
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

  unsigned int
  get_number_of_pbc() const
  {
    return number_of_declared_periodic_boundaries;
  }

  unsigned int
  get_number_of_containers() const
  {
    return boundary_combination_set_to_index.size();
  }

  std::map<std::set<types::boundary_id>, std::uint8_t>
  get_boundary_comb_to_index() const
  {
    return boundary_combination_set_to_index;
  }

  /**
 * @brief Compute the combined periodic offsets and store them in
 * combined_periodic_offsets.
 *
 * Example, we have two PBC, between boundary 0 to 1 and 6 to 3, respectively:
 *   PBC 0: (0 <-> 1)
 *   PBC 1: (6 <-> 3)
 *
 * The simulation domain is from (0.,0.) to (1., 1.) in 2D.
 * Boundary 0, X = 0
 * Boundary 1, X = 1
 * Boundary 6, Y = 0
 * Boundary 3, Y = 1
 *
 * Valid combinations generated:
 *   Set of boundary ID | Index
 *         {0}          |  0
 *         {1}          |  1
 *         {6}          |  2
 *         {3}          |  3
 *         {0,6}        |  4
 *         {0,3}        |  5
 *         {1,6}        |  6
 *         {1,3}        |  7
 *
 * Expected output:
 *  Index    |  Combined offset (Tensor)
 *    0      |          1., 0.
 *    1      |         -1., 0.
 *    2      |          0., 1.
 *    3      |          0.,-1.
 *    4      |          1., 1.
 *    5      |          1.,-1.
 *    6      |         -1., 1.
 *    7      |         -1.,-1.
 *
 */
  void
  compute_combined_periodic_offsets();

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
  set_periodic_boundaries_info(
    typename Triangulation<dim>::cell_iterator  cell,
    const unsigned int                          face_id,
    periodic_boundaries_cells_info_struct<dim> &boundaries_information);

  /**
   * @brief Check if particle is outside the cell, if so, modify the
   * location of the particle using the relevant periodic offset for the
   * boundary.
   *
   * @param[in] pb_pair_cell_info Reference to the object with periodic
   * boundary information.
   * @param[in] particles_in_pb0_cell If the particles are linked to a cell on
   * the periodic boundary 0 (true) or the periodic boundary 1 (false).
   * @param[in,out] particles_in_cell Iterator to the particles in cell.
   * @param[out] particle_has_been_moved Flag the particle has been moved to
   * the other periodic cell.
   * @param particle_to_total_periodic_displacement
   */
  void
  check_and_move_particles(
    const periodic_boundaries_cells_info_struct<dim> &pb_pair_cell_info,
    const bool                                       &particles_in_pb0_cell,
    typename Particles::ParticleHandler<dim>::particle_iterator_range
         &particles_in_cell,
    bool &particle_has_been_moved,
    typename DEM::dem_data_structures<dim>::particle_index_tensor_map
      &particle_to_total_periodic_displacement);


  /**
   * @brief Flag for periodic boundary conditions in simulation. Useful to
   * exit function when there are no periodic boundaries.
   */
  bool periodic_boundaries_enabled;



  /**
   * @brief Direction of the periodic boundaries, it is the perpendicular axis
   * of the periodic boundaries. Keys of this map are periodic boundary
   * condition indices in the .prm (subsection numbers)
   */
  std::unordered_map<unsigned int, unsigned int> directions;

  /**
   * @brief Index of the boundary conditions in the .prm (subsection numbers)
   * that correspond to periodic boundary conditions.
   */
  std::vector<unsigned int> prm_periodic_bc_index;

  /**
   * @brief A map where the keys are the boundary condition IDs of the periodic
   * boundary conditions like they are declared in the parameter file. The
   * values are the primary mesh boundary IDs associated. No need to store the
   * secondary mesh boundary IDs since the primary and secondary IDs are linked
   * on the triangulation and accessible through functions on the cells on
   * boundaries. Map key: index of BC from .prm Map value: ID of a primary
   * periodic boundary
   */
  std::unordered_map<unsigned int, types::boundary_id>
    primary_periodic_boundaries_ids;

  /**
   * @brief
   */
  std::unordered_map<unsigned int, types::boundary_id>
    secondary_periodic_boundaries_ids;

  /**
   * @brief
   */
  std::map<std::set<types::boundary_id>, std::uint8_t>
    boundary_combination_set_to_index;

  /**
   * @brief
   */
  std::vector<BoundarySet> index_to_boundary_combination;

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

  std::uint8_t number_of_declared_periodic_boundaries;
};

#endif
