// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/adaptive_sparse_contacts.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

template <int dim, typename PropertiesIndex>
AdaptiveSparseContacts<dim, PropertiesIndex>::AdaptiveSparseContacts()
  : sparse_contacts_enabled(false)
  , advect_particles_enabled(false)
{}

template <int dim, typename PropertiesIndex>
void
AdaptiveSparseContacts<dim, PropertiesIndex>::update_local_and_ghost_cell_set(
  const DoFHandler<dim> &background_dh)
{
  if (!sparse_contacts_enabled)
    return;

  local_and_ghost_cells.clear();
  for (const auto &cell : background_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        local_and_ghost_cells.insert(cell);
    }
}

template <int dim, typename PropertiesIndex>
void
AdaptiveSparseContacts<dim, PropertiesIndex>::
  calculate_granular_temperature_and_solid_fraction(
    const Particles::ParticleHandler<dim> &particle_handler,
    const std::set<typename DoFHandler<dim>::active_cell_iterator>
                   &local_and_ghost_cells_with_particles,
    Vector<double> &granular_temperature_average,
    Vector<double> &solid_fractions)
{
  // Iterating through the active cells in cell set with particles
  for (const auto &cell : local_and_ghost_cells_with_particles)
    {
      // Particles in the cell
      auto particles_in_cell = particle_handler.particles_in_cell(cell);
      const unsigned int n_particles_in_cell =
        particle_handler.n_particles_in_cell(cell);

      // Initialize variables for solid fraction computation
      double       solid_fraction = 0.0;
      double       solid_volume   = 0.0;
      const double cell_volume    = cell->measure();

      // Initialize variables for granular temperature computation
      double         granular_temperature_cell = 0.0;
      Tensor<1, dim> velocity_cell_average;
      Tensor<1, dim> cell_velocity_fluctuation_squared_average;

      // First loop over particles in cell to compute the sum of
      // particle velocity and the solid volume of the current cell
      for (auto particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          // Get particle properties
          auto particle_properties =
            particles_in_cell_iterator->get_properties();
          const double dp = particle_properties[PropertiesIndex::dp];

          for (int d = 0; d < dim; ++d)
            {
              // Get the particle velocity component (v_x, v_y & v_z if dim = 3)
              int v_axis = PropertiesIndex::v_x + d;

              // Add the velocity component value
              velocity_cell_average[d] += particle_properties[v_axis];
            }

          solid_volume += M_PI * Utilities::fixed_power<dim>(dp) / (2.0 * dim);
        }

      // Calculate average velocity in the cell (sum/n_particles)
      velocity_cell_average /= n_particles_in_cell;

      // Calculate solid fraction of cell
      solid_fraction = solid_volume / cell_volume;

      // Second loop over particle to compute the sum of the cell velocity
      // fluctuations
      for (auto particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          auto particle_properties =
            particles_in_cell_iterator->get_properties();

          for (int d = 0; d < dim; ++d)
            {
              // Get the particle velocity component (v_x, v_y & v_z if dim = 3)
              int v_axis = PropertiesIndex::v_x + d;

              cell_velocity_fluctuation_squared_average[d] +=
                Utilities::fixed_power<2>(particle_properties[v_axis] -
                                          velocity_cell_average[d]);
            }
        }

      // Calculate average granular temperature in the cell
      for (int d = 0; d < dim; ++d)
        {
          cell_velocity_fluctuation_squared_average[d] /= n_particles_in_cell;
          granular_temperature_cell +=
            cell_velocity_fluctuation_squared_average[d] / dim;
        }

      // Store the average granular temperature and solid fraction with
      // active cell index
      granular_temperature_average[cell->active_cell_index()] =
        granular_temperature_cell;
      solid_fractions[cell->active_cell_index()] = solid_fraction;
    }
}

template <int dim, typename PropertiesIndex>
void
AdaptiveSparseContacts<dim, PropertiesIndex>::identify_mobility_status(
  const DoFHandler<dim>                 &background_dh,
  const Particles::ParticleHandler<dim> &particle_handler,
  const unsigned int                     n_active_cells,
  MPI_Comm                               mpi_communicator)
{
  // If the sparse contacts are not enabled, exit the function
  if (!sparse_contacts_enabled)
    return;

  // Reset cell status containers
  cell_mobility_status.clear();

  // Reset all status to mobile if mobility status needs to be reset
  if (DEMActionManager::get_action_manager()->check_mobility_status_reset())
    {
      for (auto cell : local_and_ghost_cells)
        {
          assign_mobility_status(cell->active_cell_index(),
                                 mobility_status::mobile);
        }

      return;
    }

  // Get a copy of the active & ghost cells set to iterate over and remove cell
  // of the set when the mobility status is known. It avoids unnecessary
  // iteration for next loops. We don't want to modify the original set since
  // it is only updated when there's load balancing or reading of checkpoints.
  auto local_and_ghost_cells_copy = local_and_ghost_cells;

  // Create dummy dofs for background dof handler for the mobility_at_nodes
  const FE_Q<dim>    fe(1);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  // Get locally owned and relevant dofs
  const IndexSet &locally_owned_dofs = background_dh.locally_owned_dofs();
  IndexSet        locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(background_dh);

  // Reinit all value of mobility at nodes as inactive
  mobility_at_nodes.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           mpi_communicator);
  mobility_at_nodes = 0;

  // If the advection of particles setting is enabled (useful for CFD-DEM),
  // mobility status are different: inactive and static_active status are
  // advected and advected_active since they are not handled in the same way as
  // in DEM only. They have a special handling in the velocity integration step
  // because cell averaged velocity and acceleration are applied to particles,
  // so they need this different status. The criteria for those status are the
  // same as for DEM.
  mobility_status inactive_status = (!advect_particles_enabled) ?
                                      mobility_status::inactive :
                                      mobility_status::advected;
  mobility_status active_status   = (!advect_particles_enabled) ?
                                      mobility_status::static_active :
                                      mobility_status::advected_active;

  // Check if the cell is empty (n_particle = 0), if so, nodes and cells are
  // flagged as empty mobility status
  for (auto cell = local_and_ghost_cells_copy.begin();
       cell != local_and_ghost_cells_copy.end();)
    {
      // Check if the cell has any particles
      if (particle_handler.n_particles_in_cell(*cell) == 0)
        {
          std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
          (*cell)->get_dof_indices(local_dof_indices);

          // Assign inactive status to cell and empty status to nodes
          assign_mobility_status((*cell)->active_cell_index(),
                                 local_dof_indices,
                                 mobility_status::inactive,
                                 mobility_status::empty);

          // Remove cell from cell set and iterate to the following cell,
          // erase() function returns the next iterator
          cell = local_and_ghost_cells_copy.erase(cell);
        }
      else
        {
          // Since erase() is not called, iterator is incremented
          ++cell;
        }
    }

  // Update ghost values of mobility_at_nodes
  mobility_at_nodes.update_ghost_values();

  // Calculate the average granular temperature and solid fraction for each
  // cell currently in the local_and_ghost_cells_copy set (empty cells are
  // already removed)
  Vector<double> granular_temperature_average(n_active_cells);
  Vector<double> solid_fractions(n_active_cells);
  calculate_granular_temperature_and_solid_fraction(
    particle_handler,
    local_and_ghost_cells_copy,
    granular_temperature_average,
    solid_fractions);

  // Check if the cell is mobile by criteria:
  // * granular temperature > threshold or
  // * solid fraction of cell < threshold or
  // * is next to an empty cell
  // If so, nodes are flagged and cells are stored with mobile status
  for (auto cell = local_and_ghost_cells_copy.begin();
       cell != local_and_ghost_cells_copy.end();)
    {
      const unsigned int cell_id = (*cell)->active_cell_index();
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dof_indices);

      // Check if the cell has any empty neighbor cell by the value at nodes
      bool has_empty_node = false;
      for (auto node_id : local_dof_indices)
        {
          if (mobility_at_nodes(node_id) == mobility_status::empty)
            {
              has_empty_node = true;
              break; // No need to check other nodes
            }
        }

      // Check if the cell is mobile by criteria (granular temperature, solid
      // fraction, next to empty cell)
      if (granular_temperature_average[cell_id] >
            granular_temperature_threshold ||
          solid_fractions[cell_id] < solid_fraction_threshold || has_empty_node)
        {
          // Assign mobile status to cell and nodes
          assign_mobility_status(cell_id,
                                 local_dof_indices,
                                 mobility_status::mobile,
                                 mobility_status::mobile);

          // Remove cell from cell set
          cell = local_and_ghost_cells_copy.erase(cell);
        }
      else
        {
          ++cell;
        }
    }

  mobility_at_nodes.update_ghost_values();

  // Check if the cell is mobile by neighbor (at least a node is flagged as
  // mobile from previous check), this is the additional mobile layer.
  // If so, cells are stored in map as mobile status and nodes that are not
  // mobile are flagged as active
  for (auto cell = local_and_ghost_cells_copy.begin();
       cell != local_and_ghost_cells_copy.end();)
    {
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dof_indices);

      // Check if the cell has any mobile nodes by the value at nodes
      bool has_mobile_node = false;
      for (auto node_id : local_dof_indices)
        {
          if (mobility_at_nodes(node_id) == mobility_status::mobile)
            {
              has_mobile_node = true;
              break;
            }
        }

      if (has_mobile_node)
        {
          // Assign mobile status to cell and active to nodes
          assign_mobility_status((*cell)->active_cell_index(),
                                 local_dof_indices,
                                 mobility_status::mobile,
                                 active_status);

          // Remove cell from cell set and iterate to the following cell
          cell = local_and_ghost_cells_copy.erase(cell);
        }
      else
        {
          ++cell;
        }
    }

  mobility_at_nodes.update_ghost_values();

  // Check if the cell is active (at least a node is flagged as active from
  // previous check), this is the layer of active cells
  // If so, cells are stored with active status in the map
  for (auto cell = local_and_ghost_cells_copy.begin();
       cell != local_and_ghost_cells_copy.end();)
    {
      std::vector<types::global_dof_index> local_dofs_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dofs_indices);

      // Check if cell has active nodes
      bool has_active_nodes = false;
      for (auto node_id : local_dofs_indices)
        {
          if (mobility_at_nodes(node_id) == (int)active_status)
            {
              has_active_nodes = true;
              break;
            }
        }

      if (has_active_nodes)
        {
          // Assign active status to cell (no need to assign active to nodes)
          assign_mobility_status((*cell)->active_cell_index(), active_status);
          cell = local_and_ghost_cells_copy.erase(cell);
        }
      else
        {
          ++cell;
        }
    }

  // Store the inactive cells in the map, those are all the remaining cells
  // in the local & ghost cells set
  for (auto cell = local_and_ghost_cells_copy.begin();
       cell != local_and_ghost_cells_copy.end();
       ++cell)
    {
      // Assign inactive status to cell in map
      assign_mobility_status((*cell)->active_cell_index(), inactive_status);
    }
}

template <int dim, typename PropertiesIndex>
void
AdaptiveSparseContacts<dim, PropertiesIndex>::
  update_average_velocities_acceleration(
    Particles::ParticleHandler<dim> &particle_handler,
    const Tensor<1, 3>              &g,
    const std::vector<Tensor<1, 3>> &force,
    const double                     dt)
{
  // If the sparse contacts is disabled, exit the function.
  // Also, if the mobility status do not need to be reset, which only happens
  // at first dem iteration of the cfd iteration
  if (!sparse_contacts_enabled ||
      !DEMActionManager::get_action_manager()->check_mobility_status_reset())
    return;

  cell_velocities_accelerations.clear();

  // Tensors for velocity and acceleration * dt computation
  Tensor<1, 3> velocity_cell_average;
  Tensor<1, 3> acc_dt_cell_average;

  const Tensor<1, 3> dt_g = g * dt;

  // Loop over all the cells even if they are not mobile to reset the force
  // and torque values
  for (auto cell : local_and_ghost_cells)
    {
      if (cell->is_locally_owned())
        {
          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(cell);

          velocity_cell_average.clear();
          acc_dt_cell_average.clear();

          if (n_particles_in_cell > 0)
            {
              auto particles_in_cell = particle_handler.particles_in_cell(cell);
              for (auto &particle : particles_in_cell)
                {
                  // Get particle properties
                  auto particle_properties = particle.get_properties();
                  types::particle_index particle_id =
                    particle.get_local_index();

                  double dt_mass_inverse =
                    dt / particle_properties[PropertiesIndex::mass];

                  // Calculate the acceleration of the particle times the time
                  // step
                  Tensor<1, 3> acc_dt_particle =
                    dt_g + force[particle_id] * dt_mass_inverse;

                  // Add up the acc. * dt for the average
                  acc_dt_cell_average += acc_dt_particle;

                  for (int d = 0; d < dim; ++d)
                    {
                      // Add up the current velocity for the average velocity
                      velocity_cell_average[d] +=
                        particle_properties[PropertiesIndex::v_x + d] +
                        acc_dt_particle[d];
                    }
                }

              // Compute the average velocity and acceleration * dt
              velocity_cell_average /= n_particles_in_cell;
              acc_dt_cell_average /= n_particles_in_cell;
            }

          // Insert the average pair in map
          cell_velocities_accelerations.insert(std::make_pair(
            cell, std::make_pair(velocity_cell_average, acc_dt_cell_average)));
        }
    }
}

template class AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>;
template class AdaptiveSparseContacts<2,
                                      DEM::CFDDEMProperties::PropertiesIndex>;
template class AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>;
template class AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>;
template class AdaptiveSparseContacts<3,
                                      DEM::CFDDEMProperties::PropertiesIndex>;
template class AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>;
