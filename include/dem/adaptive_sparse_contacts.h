// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_adaptive_sparse_contacts_h
#define lethe_adaptive_sparse_contacts_h

#include <core/dem_properties.h>

#include <dem/data_containers.h>
#include <dem/dem_action_manager.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

// Special template instance for this class.
// Unsigned integer would have been a better choice, but is not working
// with the Vector class (error because of calling of abs() function)
template class LinearAlgebra::distributed::Vector<int>;

/**
 * @brief Mechanism that disables contacts computation through mobility status
 * of cells based on granular temperature.
 *
 * The general idea behind the algorithm:
 *
 * It uses the granular temperature to determinate if particles in a cell are
 * sufficiently mobile that the contact forces are worth computation.
 * A cell having a granular temperature under a threshold (default is 1e-4) will
 * have an inactive or advected status which makes them rejected in the broad
 * search step. There is also an extra layer of mobile cells around the current
 * mobile cells to propagate the motion of particles. Additionally, there is
 * another layer of cells around those mobile cells that are flagged as
 * static_active or advected_active cells.
 *
 * While the inactive status makes the particles static (no further handling),
 * the advected status will lead to computation and application of a
 * cell-averaged advected term. This term is computed from the last average
 * velocity and acceleration of the particles when the cell had a mobile status,
 * or at the first DEM time step of a CFD iteration in CFD-DEM (ensures an
 * update of the hydrodynamic forces). Those forces are reapplied after each
 * DEM time step.
 *
 * Particles in cells with an inactive status results in skipping the fine
 * contact search, the p-p & p-w forces computation and the velocity and
 * position update. It reduces the computational cost of the simulation when
 * particles are not moving.
 *
 * Particles in cells with an advected status results in skipping the fine
 * contact search, the p-p & p-w forces computations but they do have a velocity
 * and position update through the advection term. This is critical for CFD-DEM.
 *
 * In summary, cells are flagged using a mobility status:
 * - mobile (all contacts are calculated, same as if the feature is not enabled)
 * - static_active (particles with low velocity variance next to mobile
 * particles)
 * - advected_active (same as above but particles will be advected)
 * - inactive (particles next to mobile particles that allow one way contact
 * forces calculation, not computation of collision pair in the cell)
 * - advected (same as above, but particles will be advected)
 *
 * There are some edge cases that need some attention. In these cases, particles
 * in cells cannot be deactivated:
 *
 * 1. Solid fraction
 * The solid fraction of the cell is under a value (default = 40%): particles
 * in the cell are not consolidated (packed) and may move freely.
 *
 * 2. Empty cell neighbor
 * Cell having empty cell neighbor: it means that particles are next to a
 * floating wall/mesh and they may have a change of contact forces from it if
 * the wall moves or disappears.
 *
 * 3. Additional mobile cell layer
 * The cell has at least one cell neighbor flagged as mobile from the 3 criteria
 * mentioned above (granular temperature, solid fractions or floating walls):
 * this is because motion is badly propagated to particles to the cell around
 * without this additional "layer" of mobile cells.
 *
 * 4. Active cell layer
 * The cell has at least one cell neighbor flagged as mobile from the previous
 * criterion (additional layer of mobile cells): again, this is to allow the
 * motion propagation, but only the particles in contact with particles
 * from the mobile cells are considered for the contact force calculation but
 * their position is not computed at the integration step. Those cells are
 * flagged as active cells. It works with the assignment and verification of the
 * mobility status at nodes to check the status of the neighboring cells.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam PropertiesIndex Index of the properties used within the ParticleHandler.
 */
template <int dim, typename PropertiesIndex>
class AdaptiveSparseContacts
{
public:
  AdaptiveSparseContacts();

  /**
   * Mobility status flags used to identify the status at nodes and the status
   * of the cell:
   *
   * inactive (0)
   * The motion of particles in the cell is considered as negligible,
   * particles within this cell are not considered in the contact detection
   * (rejected at the broad search step), so no force calculation or
   * integration is applied
   *
   * static_active (1)
   * The motion of particles in the cell is considered as negligible, but
   * there's at least one neighbor cell that is flagged as mobile, meaning that
   * particles need to be in contact candidates lists at the broad search step,
   * particles directly in contact with the mobile cell are also considered in
   * force calculation, but none of the particles in the cell are integrated.
   *
   * advected (2)
   * The variance of the velocity of the particles in cell is low, suggesting
   * that contact forces are negligible but particles are advected by other
   * forces. Note that there are currently no mechanisms that distinct active vs
   * advected status: advection of particles is a user-defined feature. Advected
   * and active status are detected in the same way, but are not processed the
   * same at the velocity integration step. For advected particles, the last
   * computed cell-averaged velocity and acceleration is used to advect the
   * particles.
   *
   * advected_active (3)
   * Same as (1), but particles are advected by the last computed
   * cell-averaged velocity and acceleration.
   *
   * mobile (4)
   * The motion of particles in the cell is significant or there is at least
   * one cell neighbor that is mobile by criteria (see the
   * identify_mobility_status() description), particles need to be in contact
   * candidate lists at the broad search step and particles are treated as
   * usual (force calculation and integration)
   *
   * empty (5)
   * This status is only used for the node-based mobility status identification,
   * no cells are flagged as empty, only node can be identified as empty.
   * Without this identification of the empty cells, we cannot identify the cell
   * that have an empty neighbor cell, which is critical for simulations using
   * floating walls or mesh.
   */

  enum mobility_status : unsigned int
  {
    inactive        = 0, // used for cells
    static_active   = 1, // used for cells and nodes
    advected        = 2, // used for cells and nodes
    advected_active = 3, // used for cells and nodes
    mobile          = 4, // used for cells and nodes
    empty           = 5  // used for nodes only
  };


  /**
   * @brief Set the threshold values for the mobile status criteria (granular
   * temperature and solid fraction) and the flag for the advection of
   * particles.
   *
   * @param[in] granular_temperature The threshold value for the granular
   * temperature.
   * @param[in] solid_fraction The threshold value for the volumic solid
   * fraction.
   * @param[in] advect_particles The flag for the advection of particles.
   */
  inline void
  set_parameters(const double granular_temperature,
                 const double solid_fraction,
                 const bool   advect_particles)
  {
    // If the function is reached, the adaptive sparse contacts is enabled.
    sparse_contacts_enabled = true;

    // Communicate to the action manager that the sparse contacts is enabled
    DEMActionManager::get_action_manager()->set_sparse_contacts_enabled();

    granular_temperature_threshold = granular_temperature;
    solid_fraction_threshold       = solid_fraction;
    advect_particles_enabled       = advect_particles;
  }

  /**
   * @brief Create or update a set of the active and ghost cells so that
   * there is no loop over all the cells of the triangulation for the granular
   * temperature and solid fraction calculation, and during the identification
   * of the mobility status. This set prevents 4 iteration steps over all the
   * cells and the verification if the cell is locally owned, ghost or not.
   * This set is updated at every load balance step since cells are
   * redistributed among processors.
   *
   * @param[in] background_dh The DoFHandler of the background grid.
   */
  void
  update_local_and_ghost_cell_set(const DoFHandler<dim> &background_dh);

  /**
   * @brief Calculate or update the cell-averaged velocities and accelerations
   * of the cells. This is used for the advection of particles during
   * integration and the velocities and accelerations map is fully updated at
   * least at each CFD time step since the first DEM iteration computes every
   * contacts. This map is also partially updated when cell is mobile in the
   * Verlet integration (not implemented for other scheme).
   * Since this step loops over all the particles of mobile cells and new values
   * of velocities are computed, this helps to have a more accurate value of the
   * average velocities and accelerations when applied to the particles in
   * next DEM time steps.
   *
   * @param[in] particle_handler The particle handler that contains all the
   * particles.
   * @param[in] g The external forces, same for all cells (gravity or other
   * sources).
   * @param[in] force The contact or hydrodynamic forces.
   * @param[in] dt The DEM time step.
   */
  void
  update_average_velocities_acceleration(
    Particles::ParticleHandler<dim> &particle_handler,
    const Tensor<1, 3>              &g,
    const std::vector<Tensor<1, 3>> &force,
    const double                     dt);

  /**
   * @brief Identify of the mobility status of each cell through a node-based
   * identification and check. Only the active and ghost cells are processed.
   *
   * The following 4 checks (search loops) are done:
   *
   * 1. Check if the cell is empty (n_particle = 0), if so, nodes and cells are
   * flagged as empty mobility status (5)
   *
   * 2. Check if the cell is mobile by criteria (average granular temperature >
   * threshold, solid fraction < threshold or has at least one empty node from
   * previous check), if so, nodes are flagged and cells are stored with mobile
   * mobility status (4)
   *
   * 3. Check if the cell is mobile by neighbor (at least a node is flagged as
   * mobile from previous check), if so, cells are stored in map as mobile
   * status (4) and nodes that are not mobile are flagged as active (1/3)
   *
   * 4. Check if the cell is active (at least a node is flagged as active from
   * previous check), if so, cells are stored with active status in the map
   * (1/3)
   *
   * The remaining cells are inactive (0)
   *
   * @param[in] background_dh The dof handler of the background grid.
   * @param[in] particle_handler The particle handler that contains all the
   * particles.
   * @param[in] n_active_cells The number of active cells in triangulation.
   * @param[in] mpi_communicator The MPI communicator.
   */
  void
  identify_mobility_status(
    const DoFHandler<dim>                 &background_dh,
    const Particles::ParticleHandler<dim> &particle_handler,
    const unsigned int                     n_active_cells,
    MPI_Comm                               mpi_communicator);

  /**
   * @brief Map the periodic nodes pairs of the triangulation using the
   * constraints. It allows to compare the mobility status of the nodes with the
   * status of the periodic node.
   *
   * Note: this might not be the efficient way to pair the periodic nodes.
   *
   * @param constraints[in] The background constraints of triangulation.
   */
  inline void
  map_periodic_nodes(const AffineConstraints<double> &constraints)
  {
    periodic_node_ids.clear();

    IndexSet local_lines = constraints.get_local_lines();
    for (auto i : local_lines)
      {
        for (auto j : local_lines)
          {
            // Map the values of the periodic nodes (with repetition)
            if (constraints.are_identity_constrained(i, j))
              {
                periodic_node_ids.insert(std::make_pair(i, j));
              }
          }
      }
  }

  /**
   * @brief Find the mobility status of a cell.
   *
   * @param[in] cell The iterator of the cell that needs mobility evaluation.
   */
  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    return cell_mobility_status.at(cell->active_cell_index());
  }

  /**
   * @brief Convert the map of mobility status to a vector of mobility status
   * because map cannot be used as is in the pvd post-processing or any data
   * out, it needs to be converted to a vector of mobility status by active cell
   * index.
   *
   * @param[out] status The initiated vector for the conversion.
   */
  void
  get_mobility_status_vector(Vector<float> &status)
  {
    for (auto &cell_to_status : cell_mobility_status)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  }

  /**
   * @brief Give the map of the mobility status of the cells.
   */
  typename DEM::dem_data_structures<dim>::cell_index_int_map &
  get_mobility_status()
  {
    return cell_mobility_status;
  }

  /**
   * @brief Give the map of the cell-averaged velocities and accelerations * dt.
   */
  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::pair<Tensor<1, 3>, Tensor<1, 3>>> &
  get_velocities_accelerations()
  {
    return cell_velocities_accelerations;
  }

  /**
   * @brief Give the advected particles flag.
   */
  bool
  has_advected_particles() const
  {
    return advect_particles_enabled;
  }

private:
  /**
   * @brief Calculate the granular temperature and solid fraction approximation
   * (pcm method). Those values are criteria for cell mobility.
   *
   * @param[in] particle_handler The particle handler that contains all the
   * particles.
   * @param[in] local_and_ghost_cells_with_particles The set of locally owned
   * and ghost cells that empty cells were removed from.
   * @param[out] cell_granular_temperature The empty vector of granular
   * temperature.
   * @param[out] cell_solid_fraction The empty vector of solid fraction.
   */
  void
  calculate_granular_temperature_and_solid_fraction(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::set<typename DoFHandler<dim>::active_cell_iterator>
                   &local_and_ghost_cells_with_particles,
    Vector<double> &cell_granular_temperature,
    Vector<double> &cell_solid_fraction);

  /**
   * @brief Assign the mobility status to the cell, its nodes, and its periodic
   * coinciding nodes from the periodic nodes map. If required, this will
   * prevent overwriting a node status prevailing over the initial one.
   * For instance, empty status must not be overwritten by mobile status because
   * all cells with this node will be part of the additional mobile layer and
   * will wrongly propagate mobility status to the next cells.
   *
   * @param[in] cell_id The current cell index.
   * @param[in] dof_indices The vector of the DoF indices of the cell.
   * @param[in] cell_status The mobility status of cell to assign.
   * @param[in] node_status The current mobility status of the node.
   */
  inline void
  assign_mobility_status(unsigned int                          cell_id,
                         std::vector<types::global_dof_index> &dof_indices,
                         const int                             cell_status,
                         const int                             node_status)
  {
    cell_mobility_status.insert({cell_id, cell_status});

    // Assign mobility status to nodes without overwriting empty or mobile nodes
    // in regard to the case
    for (auto node_id : dof_indices)
      {
        // Prevailing mobility status of the node and assignation
        int status_assigned = std::max(node_status, mobility_at_nodes(node_id));
        mobility_at_nodes(node_id) = status_assigned;

        // Check if node has a periodic node and assign the same mobility status
        // if prevailing on the one on the periodic node
        auto it = periodic_node_ids.find(node_id);
        if (it != periodic_node_ids.end())
          {
            mobility_at_nodes(it->second) =
              std::max(status_assigned, mobility_at_nodes(it->second));
          }
      }
  }

  /**
   * @brief Assign the mobility status to the cell. No need to update the node
   * status. Only use for active status assignment since nodes are not evaluated
   * after this step.
   *
   * @param[in] cell_id The current cell index.
   * @param[in] cell_status The mobility status of cell to assign.
   */
  inline void
  assign_mobility_status(unsigned int cell_id, const int cell_status)
  {
    cell_mobility_status.insert({cell_id, cell_status});
  }

  /**
   * @brief Set of locally owned and ghost cells: <local/ghost cells>
   * Used to loop over only the locally owned and ghost cells without looping
   * over all the cells in the triangulation numerous times.
   */
  std::set<typename DoFHandler<dim>::active_cell_iterator>
    local_and_ghost_cells;

  /**
   * @brief Map of cell mobility status: <cell index: mobility status>
   */
  typename DEM::dem_data_structures<dim>::cell_index_int_map
    cell_mobility_status;

  /**
   * @brief Vector of mobility status at nodes: [mobility status]
   * Used to check the value at node to determine the mobility status of the
   * cell, this type of vector is used to allow update values in parallel.
   */
  LinearAlgebra::distributed::Vector<int> mobility_at_nodes;

  /**
   * @brief Map of periodic nodes: <periodic node index: coinciding node index>
   */
  std::unordered_map<unsigned int, unsigned int> periodic_node_ids;

  /**
   * @brief Map of cell velocities and accelerations * dt:
   * <cell iterator: <velocity, acceleration * dt>>
   */
  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::pair<Tensor<1, 3>, Tensor<1, 3>>>
    cell_velocities_accelerations;

  /**
   * @brief Flag for adaptive sparse contacts enabling. Useful to exit functions.
   */
  bool sparse_contacts_enabled;

  /**
   * @brief Flag for the advection of particles.
   */
  bool advect_particles_enabled;

  /**
   * @brief Threshold value for granular temperature.
   */
  double granular_temperature_threshold;

  /**
   * @brief Threshold value for solid fraction.
   */
  double solid_fraction_threshold;
};

#endif // lethe_adaptive_sparse_contacts_h
