/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2022 by the Lethe authors
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
 */

#include <core/dem_properties.h>

#include <dem/data_containers.h>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_vector.templates.h>

#include <deal.II/particles/particle_handler.h>

#include <vector>

using namespace dealii;

#ifndef lethe_disable_contacts_h
#  define lethe_disable_contacts_h

// Special template instance for this class.
// Unsigned integer would have been a better choice, but is not working
// with the Vector class (error because of calling of abs() function)
template class LinearAlgebra::distributed::Vector<int>;

/**
 * The general idea behind the algorithm:
 * It uses the granular temperature to determinate if particles in a cell are
 * mobile enough that the contact forces are worth computation. A cell having a
 * mobile enough that the contact forces are worth computation. A cell having a
 * granular temperature under a value (default is 1e-4) will have an inactive
 * status which makes them rejected in the broad search step. This results in no
 * contact forces or velocity integration in future steps for those particles,
 * which make less computation cost for the simulation.
 *
 * Cells may be flagged as that called mobility status:
 * - mobile (everything is calculated as if feature is not enabled)
 * - active (particles with low motion next to mobile particles)
 * - inactive (particles are not in a neighbor’s candidate list of cells around)
 *
 * There are some edge cases that need some attention. In these cases, particles
 * in cells can't be deactivated:
 *
 * 1. The solid fraction of the cell is under a value (default = 40%): particles
 * in the cell may have other forces that can move them and it is supposed that
 * with this fraction, there is not enough particles around to hold them in
 * their position.
 *
 * 2. Cells having empty cell neighbors: it means that particles are next to a
 * floating wall/mesh and they may have a change of contact forces from it if
 * the wall disappears or is moving.
 *
 * 3. The cell has cell neighbors which is flagged as mobile from the 3 criteria
 * mentioned above (granular temperature, solid fractions or floating walls):
 * this is because motion is badly propagated to particles to the cell around
 * without this additional "layer" of mobile cells.
 *
 * 4. The cell has cell neighbors which is flagged as mobile from the previous
 * criterion (additional layer of mobile cells): this is again to allow the
 * propagation of motion, but only the particles in contact with particles
 * from the mobile cells are considered for the contact force calculation but
 * their position is not computed at the integration step. Those cells are
 * flagged as active cells. It works with the assignment and verification of the
 * mobility status at nodes to check the status of the neighboring cells.
 *
 */
template <int dim>
class DisableContacts
{
public:
  DisableContacts<dim>();

  /**
   * Mobility status flag used to identify the status at nodes and the status
   * of the cell:
   *
   * inactive (0)
   * The movement of particles in the cell is considered as negligible,
   * particles within this cell are not considered in the contact detection
   * (rejected at the broad search step), so no force calculation or
   * integration is applied
   *
   * active (1)
   * The movement of particles in the cell is considered as negligible, but
   * there's at least one neighbor cell that is flagged as mobile, meaning that
   * particles need to be in contact candidates lists at the broad search step,
   * particles directly in contact with the mobile cell are also considered in
   * force calculation, but none of the particles in the cell are integrated.
   *
   * mobile (4)
   * The movement of particles in the cell is significant or there is at least
   * one neighbor cell that is mobile by criteria (see the
   * identify_mobility_status() description), particles need to be in contact
   * candidates lists at the broad search step and particles are treated as
   * usual (force calculation and integration)
   *
   * empty (5)
   * This status is only used for the node-based mobility status identification,
   * no cells are flagged as empty, only node can by identify as empty. Without
   * this identification of the empty cells, we can't identify the cell that
   * have a empty neighbor cell, which is critical for simulations using
   * floating walls or mesh
   */

  // TODO : add doc for advected and advected_active
  enum mobility_status : unsigned int
  {
    inactive        = 0, // used for cells
    active          = 1, // used for cells and nodes
    advected        = 2, // used for cells and nodes
    advected_active = 3, // used for cells and nodes
    mobile          = 4, // used for cells and nodes
    empty           = 5  // used for nodes only
  };

  /**
   * @brief Create or update a set of the active and ghost cells so that we don't
   * have to loop over all the cells in the triangulation for the granular
   * temperature and solid fraction calculation, and during the identification
   * of the mobility status. This set prevent 4 iteration steps over all the
   * cells  + the verification if the cell is locally owned, ghost or not.
   * This set is updated at every load balance step since cells are
   * redistributed among processors.
   *
   * @param background_dh The DoFHandler of the background grid
   */
  void
  update_local_and_ghost_cell_set(const DoFHandler<dim> &background_dh);

  /**
   * @brief Carries out the identification of the mobility status of each cell
   * through a node-based identification and check. Only the active and ghost
   * cells are processed.
   *
   * The following 4 checks (search loops) are done:
   *
   * 1. Check if the cell is empty (n_particle = 0), if so, nodes and cells are
   * flagged as empty mobility status (3)
   *
   * 2. Check if the cell is mobile by criteria (average granular temperature >
   * threshold, solid fraction < threshold or has at least one empty node from
   * previous check), if so, nodes are flagged and cells are stored with mobile
   * mobility status (2)
   *
   * 3. Check if the cell is mobile by neighbor (at least a node is flagged as
   * mobile from previous check), if so, cells are stored in map as mobile
   * status (2) and nodes that are not mobile are flagged as active (1)
   *
   * 4. Check if the cell is active (at least a node is flagged as active from
   * previous check), if so, cells are stored with active status in the map (1)
   *
   * The remaining cells are inactive (0)
   *
   * @param background_dh The dof handler of the background grid
   *
   * @param particle_handler The particle handler that contains all the particles
   *
   * @param mpi_communicator The MPI communicator
   */
  void
  identify_mobility_status(
    const DoFHandler<dim> &                background_dh,
    const Particles::ParticleHandler<dim> &particle_handler,
    const unsigned int                     n_active_cells,
    MPI_Comm                               mpi_communicator);

  inline void
  map_periodic_nodes(const AffineConstraints<double> &constraints)
  {
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

  inline void
  assign_to_periodic_node(const unsigned int node_id,
                          const unsigned int mobility_status,
                          const bool         check_max_mobility_status = false)
  {
    // Check if node has a periodic node
    auto it = periodic_node_ids.find(node_id);
    if (it != periodic_node_ids.end())
      {
        if (!check_max_mobility_status)
          {
            // Assign the mobility status to the periodic node
            mobility_at_nodes(it->second) = mobility_status;
          }
        else
          {
            // Check if the mobility status is higher than the one already
            // assigned to the periodic node (prevents from overwriting)
            mobility_at_nodes(it->second) =
              std::max((int)mobility_status, mobility_at_nodes(it->second));
          }
      }
  }


  /**
   * @brief Find the mobility status of a cell
   *
   * @param cell The iterator of the cell that needs mobility evaluation
   */
  inline unsigned int
  check_cell_mobility(
    const typename Triangulation<dim>::active_cell_iterator &cell) const
  {
    return cell_mobility_status.at(cell->active_cell_index());
  }

  /**
   * @brief Set the threshold values for the mobile status criteria (granular
   * temperature and solid fraction)
   *
   * @param granular_temperature The threshold value for the granular temperature
   *
   * @param solid_fraction The threshold value for the solid fraction (volume of
   * particles in the cell)
   */
  void
  set_threshold_values(const double granular_temperature,
                       const double solid_fraction,
                       const double advect_particles)
  {
    granular_temperature_threshold = granular_temperature;
    solid_fraction_threshold       = solid_fraction;
    advect_particles_enabled       = advect_particles;
  }

  /**
   * @brief Convert the map of mobility status to a vector of mobility status
   * because map can't be used as is in the pvd post-processing or any data out,
   * it needs to be converted to a vector of mobility status by active cell
   * index
   *
   * @param status The initiated vector for the conversion
   */
  void
  get_mobility_status_vector(Vector<float> &status)
  {
    for (auto &cell_to_status : cell_mobility_status)
      {
        status[cell_to_status.first] = cell_to_status.second;
      }
  };


  LinearAlgebra::distributed::Vector<int>
  get_mobility_output()
  {
    return mobility_at_nodes;
  }

  typename DEM::dem_data_structures<dim>::cell_index_int_map &
  get_mobility_status()
  {
    return cell_mobility_status;
  }

  void
  update_average_velocities_acceleration(
    Particles::ParticleHandler<dim> &particle_handler,
    const Tensor<1, 3> &             g,
    const std::vector<Tensor<1, 3>> &force,
    const double                     dt)
  {
    for (auto &cell : this->local_and_ghost_cells)
      {
        // We loop over the particles, even if cell is not mobile, to reset
        // the force and torques value of the particle with the particle id
        auto particles_in_cell = particle_handler.particles_in_cell(cell);
        const unsigned int n_particles_in_cell =
          particle_handler.n_particles_in_cell(cell);

        if (n_particles_in_cell > 0)
          {
            Tensor<1, 3> velocity_cell_average;
            Tensor<1, 3> acceleration_dt;

            for (auto &particle : particles_in_cell)
              {
                // Get particle properties
                auto &particle_properties         = particle.get_properties();
                types::particle_index particle_id = particle.get_local_index();

                for (int d = 0; d < dim; ++d)
                  {
                    // Get the particle velocity component (v_x, v_y & v_z if
                    // dim = 3)
                    int v_axis = DEM::PropertiesIndex::v_x + d;
                    velocity_cell_average[d] += particle_properties[v_axis];
                  }

                acceleration_dt +=
                  force[particle_id] /
                    particle_properties[DEM::PropertiesIndex::mass] +
                  g;
              }

            velocity_cell_average /= n_particles_in_cell;
            acceleration_dt *= dt / n_particles_in_cell;

            // Update acceleration for the mobile cell only
            cell_velocities_accelerations[cell] = {velocity_cell_average,
                                                   acceleration_dt};
          }
      }
  }

  std::vector<Tensor<1, 3>>
  get_cell_velocities(unsigned int n_active_cell)
  {
    std::vector<Tensor<1, 3>> cell_velocities(n_active_cell);

    for (auto &cell : this->local_and_ghost_cells)
      {
        cell_velocities[cell->active_cell_index()] = cell_velocities_map[cell];
      }

    return cell_velocities;
  }

  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::pair<Tensor<1, 3>, Tensor<1, 3>>> &
  get_velocities_accelerations()
  {
    return cell_velocities_accelerations;
  }


  void
  calculate_cell_acceleration(Particles::ParticleHandler<dim> &particle_handler,
                              const Tensor<1, 3> &             g,
                              std::vector<Tensor<1, 3>> &      force,
                              bool only_mobile = false)
  {
    for (auto &cell : this->local_and_ghost_cells)
      {
        Tensor<1, 3> empty_cell_acceleration({0.0, 0.0, 0.0});
        if (only_mobile)
          if (cell_mobility_status[cell->active_cell_index()] !=
              mobility_status::mobile)
            continue;

        // We loop over the particles, even if cell is not mobile, to reset
        // the force and torques value of the particle with the particle id
        auto particles_in_cell = particle_handler.particles_in_cell(cell);
        const unsigned int n_particles_in_cell =
          particle_handler.n_particles_in_cell(cell);

        if (n_particles_in_cell > 0)
          {
            Tensor<1, 3> acceleration;

            for (auto &particle : particles_in_cell)
              {
                // Get particle properties
                auto &particle_properties         = particle.get_properties();
                types::particle_index particle_id = particle.get_local_index();
                acceleration +=
                  force[particle_id] /
                    particle_properties[DEM::PropertiesIndex::mass] +
                  g;
              }

            acceleration /= n_particles_in_cell;

            // Update acceleration for the mobile cell only
            cell_acceleration_map[cell] = acceleration;
          }
      }
  }

  void get_cell_acceleration(std::vector<Tensor<1, 3>> &acceleration)
  {
    for (auto &cell : this->local_and_ghost_cells)
      {
        acceleration[cell->active_cell_index()] = cell_acceleration_map[cell];
      }
  }

  void
  get_acceleration_vector(Vector<float> &acceleration_x,
                          Vector<float> &acceleration_y,
                          Vector<float> &acceleration_z)
  {
    for (auto &cell : this->local_and_ghost_cells)
      {
        acceleration_x[cell->active_cell_index()] =
          cell_acceleration_map[cell][0];
        acceleration_y[cell->active_cell_index()] =
          cell_acceleration_map[cell][1];
        acceleration_z[cell->active_cell_index()] =
          cell_acceleration_map[cell][2];
      }
  };

  bool
  has_advected_particles() const
  {
    return advect_particles_enabled;
  }

private:
  /**
   * @brief Carries out the calculation of the granular temperature and solid
   * fraction approximation (pcm method) in each active cell. Those values are
   * criteria for cell mobility
   *
   * solid fraction = n particles * particle volume / cell volume
   *
   * Granular temperature:
   * 1. average velocity of cell = sum of particle velocity /
   * n_particles
   * 2. average of cell velocity fluctuation squared = sum of particle velocity
   * fluctuation squared / n_particles
   * 3. granular temperature = sum of average of cell velocity fluctuation
   * squared / dim
   *
   * @param particle_handler The particle handler that contains all the particles
   *
   * @param current_local_and_ghost_cells The set of locally owned and ghost
   * cells that empty cells are removed from
   *
   * @param cell_granular_temperature The empty vector of granular temperature
   *
   * @param cell_solid_fraction The empty vector of solid fraction
   */
  void
  calculate_granular_temperature_and_solid_fraction(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::set<typename DoFHandler<dim>::active_cell_iterator>
      &             local_and_ghost_cells_with_particles,
    Vector<double> &cell_granular_temperature,
    Vector<double> &cell_solid_fraction);

  // Map of cell mobility status, the key is the active cell index and the value
  // is the mobility status
  typename DEM::dem_data_structures<dim>::cell_index_int_map
    cell_mobility_status;

  // Set of locally owned and ghost cells, used to loop over only the locally
  // owned and ghost cells without looping over all the cells in the
  // triangulation many times
  std::set<typename DoFHandler<dim>::active_cell_iterator>
    local_and_ghost_cells;

  TrilinosWrappers::MPI::Vector mobility_output;

  // Vector of mobility status at nodes, used to check the value at node to
  // determine the mobility status of the cell, this type of vector is used
  // to allow update values in parallel
  LinearAlgebra::distributed::Vector<int> mobility_at_nodes;

  std::unordered_map<unsigned int, unsigned int> periodic_node_ids;

  // Particle advection
  bool advect_particles_enabled;

  // Threshold values for granular temperature and solid fraction
  double granular_temperature_threshold;
  double solid_fraction_threshold;


  std::vector<Tensor<1, 3>> cell_acceleration;

  std::map<typename DoFHandler<dim>::active_cell_iterator, Tensor<1, 3>>
    cell_acceleration_map;

  std::map<typename DoFHandler<dim>::active_cell_iterator, Tensor<1, 3>>
    cell_velocities_map;

  std::map<typename Triangulation<dim>::active_cell_iterator,
           std::pair<Tensor<1, 3>, Tensor<1, 3>>>
    cell_velocities_accelerations;
};

#endif // lethe_disable_contacts_h
