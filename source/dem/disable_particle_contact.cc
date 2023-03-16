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

#include <dem/disable_particle_contact.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

template <int dim>
DisableParticleContact<dim>::DisableParticleContact()
{}

template <int dim>
void
DisableParticleContact<dim>::update_active_ghost_cell_set(
  const DoFHandler<dim> &background_dh)
{
  active_ghost_cells.clear();
  for (const auto &cell : background_dh.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          active_ghost_cells.insert(cell);
        }
    }
}

template <int dim>
void
DisableParticleContact<dim>::calculate_granular_temperature_solid_fraction(
  const Particles::ParticleHandler<dim> &particle_handler,
  Vector<double> &                       granular_temperature_average,
  Vector<double> &                       solid_fractions)
{
  // Iterating through the active cells in the active and ghost cells set
  for (const auto &cell : active_ghost_cells)
    {
      double granular_temperature_cell = 0.0;
      double solid_fraction            = 0.0;

      // Particles in the cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
                         particles_in_cell = particle_handler.particles_in_cell(cell);
      const unsigned int n_particles_in_cell =
        particle_handler.n_particles_in_cell(cell);

      // Check if the cell has any particles
      if (n_particles_in_cell > 0)
        {
          // Initialize variables for average velocity
          Tensor<1, dim> velocity_cell_sum;
          Tensor<1, dim> velocity_cell_average;

          // Initialize variables for solid fraction
          double       solid_volume = 0.0;
          const double cell_volume  = cell->measure();

          // Initialize velocity fluctuations
          Tensor<1, dim> cell_velocity_fluctuation_squared_sum;
          Tensor<1, dim> cell_velocity_fluctuation_squared_average;

          // First loop over particles in cell to calculation the average
          // velocity and the void fraction
          for (typename Particles::ParticleHandler<dim>::
                 particle_iterator_range::iterator particles_in_cell_iterator =
                   particles_in_cell.begin();
               particles_in_cell_iterator != particles_in_cell.end();
               ++particles_in_cell_iterator)
            {
              auto &particle_properties =
                particles_in_cell_iterator->get_properties();

              for (int d = 0; d < dim; ++d)
                {
                  velocity_cell_sum[d] +=
                    particle_properties[DEM::PropertiesIndex::v_x + d];
                }

              solid_volume +=
                M_PI * pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
                (2.0 * dim);
            }

          // Calculate average velocity in the cell
          for (int d = 0; d < dim; ++d)
            velocity_cell_average[d] =
              velocity_cell_sum[d] / n_particles_in_cell;

          // Calculate solid fraction of cell
          solid_fraction = solid_volume / cell_volume;

          // Second loop over particle to calculate the average granular
          // temperature
          for (typename Particles::ParticleHandler<dim>::
                 particle_iterator_range::iterator particles_in_cell_iterator =
                   particles_in_cell.begin();
               particles_in_cell_iterator != particles_in_cell.end();
               ++particles_in_cell_iterator)
            {
              auto &particle_properties =
                particles_in_cell_iterator->get_properties();

              for (int d = 0; d < dim; ++d)
                {
                  cell_velocity_fluctuation_squared_sum[d] +=
                    (particle_properties[DEM::PropertiesIndex::v_x + d] -
                     velocity_cell_average[d]) *
                    (particle_properties[DEM::PropertiesIndex::v_x + d] -
                     velocity_cell_average[d]);
                }
            }

          // Calculate average granular temperature in the cell
          for (int d = 0; d < dim; ++d)
            {
              cell_velocity_fluctuation_squared_average[d] =
                cell_velocity_fluctuation_squared_sum[d] / n_particles_in_cell;
              granular_temperature_cell +=
                (1.0 / dim) * cell_velocity_fluctuation_squared_average[d];
            }
        }

      // Store the average granular temperature and solid fraction with
      // active cell index
      granular_temperature_average[cell->active_cell_index()] =
        granular_temperature_cell;
      solid_fractions[cell->active_cell_index()] = solid_fraction;
    }
}

template <int dim>
void
DisableParticleContact<dim>::identify_mobility_status(
  const DoFHandler<dim> &                background_dh,
  const Particles::ParticleHandler<dim> &particle_handler,
  const unsigned int                     n_active_cells,
  MPI_Comm                               mpi_communicator)
{
  // Reset cell status containers
  cell_mobility_status_map.clear();

  // Get a copy of the active & ghost cells set to iterate over and remove cell
  // of the set when the mobility status is known to avoid unnecessary
  // iterations for next loops. We don't want to modify the original set since
  // it is only updated when there's load balancing or reading of checkpoints.
  auto active_ghost_cells_copy = active_ghost_cells;

  // Create dummy dofs for background dof handler for the mobility_at_nodes
  const FE_Q<dim>    fe(1);
  const unsigned int dofs_per_cell = fe.dofs_per_cell;

  // Get locally owned and relevant dofs
  const IndexSet locally_owned_dofs = background_dh.locally_owned_dofs();
  IndexSet       locally_relevant_dofs =
    DoFTools::extract_locally_relevant_dofs(background_dh);

  // Reinit all value of mobility at nodes as inactive (0)
  mobility_at_nodes.reinit(locally_owned_dofs,
                           locally_relevant_dofs,
                           mpi_communicator);
  mobility_at_nodes = 0;

  // Calculate the average granular temperature and solid fraction for each
  // cells
  Vector<double> granular_temperature_average(n_active_cells);
  Vector<double> solid_fractions(n_active_cells);
  calculate_granular_temperature_solid_fraction(particle_handler,
                                                granular_temperature_average,
                                                solid_fractions);

  // Check if the cell is empty (n_particle = 0), if so, nodes and cells are
  // flagged as empty mobility status (3) (empty cells are not stored in the
  // map).
  for (auto cell = active_ghost_cells_copy.begin();
       cell != active_ghost_cells_copy.end();)
    {
      // Check if the cell has any particles
      if (particle_handler.n_particles_in_cell(*cell) == 0)
        {
          std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
          (*cell)->get_dof_indices(local_dof_indices);

          // Remove cell from cell set and iterate to the following cell, the
          // erase function returns the next iterator.
          cell = active_ghost_cells_copy.erase(cell);

          // Assign empty status to all nodes
          for (auto node_id : local_dof_indices)
            {
              mobility_at_nodes(node_id) = mobility_status::empty;
            }
        }
      else
        {
          // Since erase() is not called, we need to increment the iterator
          ++cell;
        }
    }

  // Update ghost values of mobility_at_nodes
  mobility_at_nodes.update_ghost_values();

  // Check if the cell is mobile by criteria:
  // * granular temperature > threshold or
  // * solid fraction of cell < threshold or
  // * is next to an empty cell
  // If so, nodes are flagged and cells are stored with mobile status (2)
  for (auto cell = active_ghost_cells_copy.begin();
       cell != active_ghost_cells_copy.end();)
    {
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dof_indices);

      // Check if the cell has any empty neighbor cell by the value at nodes
      bool has_empty_neighbor = false;
      for (auto node_id : local_dof_indices)
        {
          if (mobility_at_nodes[node_id] == mobility_status::empty)
            {
              has_empty_neighbor = true;
              break; // No need to check the other nodes
            }
        }

      // Check if the cell is mobile by criteria
      // (granular temperature or solid fraction or empty neighbor)
      const unsigned int cell_id = (*cell)->active_cell_index();
      if (granular_temperature_average[cell_id] >
            granular_temperature_threshold ||
          solid_fractions[cell_id] < solid_fraction_threshold ||
          has_empty_neighbor)
        {
          // Assign mobile status to cell in map
          cell_mobility_status_map.insert({cell_id, mobility_status::mobile});

          // Remove cell from cell set and iterate to the following cell
          cell = active_ghost_cells_copy.erase(cell);

          // Assign mobile status to nodes but don't overwrite empty nodes.
          // This prevents assigning mobile status to the empty cells in the
          // next check (additional mobile layer) since we don't verify again if
          // cell is empty or not.
          for (auto node_id : local_dof_indices)
            {
              // Possible cases: (nodes are initialized to inactive (0))
              // node = max(mobile (2), inactive (0)) = mobile (2)
              // node = max(mobile (2), empty (3))    = empty (3)
              mobility_at_nodes(node_id) =
                std::max((int)mobility_status::mobile,
                         mobility_at_nodes[node_id]);
            }
        }
      else
        {
          // Since erase() is not called, we need to increment the iterator
          ++cell;
        }
    }

  mobility_at_nodes.update_ghost_values();

  // Check if the cell is mobile by neighbor (at least a node is flagged as
  // mobile from previous check), this is the additional mobile layer.
  // If so, cells are stored in map as mobile status (2) and nodes that are not
  // mobile are flagged as active (1)
  for (auto cell = active_ghost_cells_copy.begin();
       cell != active_ghost_cells_copy.end();)
    {
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dof_indices);

      // Check if the cell has any mobile nodes to know if we need to iterate
      // over the next cell (no call of erase() function)
      bool has_mobile_node = false;

      for (auto node_id : local_dof_indices)
        {
          // Check if node is mobile and assign mobile status to the cell
          if (mobility_at_nodes[node_id] == mobility_status::mobile)
            {
              // Assign mobile status to cell in map
              const unsigned int cell_id = (*cell)->active_cell_index();
              cell_mobility_status_map.insert(
                {cell_id, mobility_status::mobile});

              // Remove cell from cell set and iterate to the following cell,
              // also label it as mobile to avoid double iteration
              cell            = active_ghost_cells_copy.erase(cell);
              has_mobile_node = true;

              // Assign active status to nodes except mobile because
              // this will cause to propagate the mobile status to the
              // neighbors in this loop since the mobility check at node
              // is executed in the same container that we are assigning new
              // mobility status
              for (auto node_id : local_dof_indices)
                {
                  // Possible cases: (nodes are initialized to inactive (0))
                  // node = max(active (1), inactive (0)) = active (1)
                  // node = max(active (1), empty (3))    = empty (3)
                  // node = max(active (1), mobile (2))   = mobile (2)
                  mobility_at_nodes[node_id] =
                    std::max((int)mobility_status::active,
                             mobility_at_nodes[node_id]);
                }
              break; // No need to check the other nodes
            }
        }

      // Since erase() is not called (no mobile node), we need to increment the
      // iterator
      if (!has_mobile_node)
        {
          ++cell;
        }
    }

  mobility_at_nodes.update_ghost_values();

  // Check if the cell is active (at least a node is flagged as active from
  // previous check), this is the layer of active cells
  // If so, cells are stored with active status in the map (1)
  for (auto cell = active_ghost_cells_copy.begin();
       cell != active_ghost_cells_copy.end();)
    {
      std::vector<types::global_dof_index> local_dofs_indices(dofs_per_cell);
      (*cell)->get_dof_indices(local_dofs_indices);


      bool has_active_nodes = false;
      bool has_mobile_nodes = false;

      // Check if cell has active and/or mobile nodes
      for (auto node_id : local_dofs_indices)
        {
          has_active_nodes =
            (mobility_at_nodes[node_id] == mobility_status::active) ||
            has_active_nodes;

          has_mobile_nodes =
            (mobility_at_nodes[node_id] == mobility_status::mobile) ||
            has_mobile_nodes;
        }

      // Active nodes with mobile nodes means that the cell is part of the
      // additional mobile layer.
      // Active nodes and no mobile nodes means that the cell is a neighbor
      // of the additional layer of mobile cells, this is an active cell.
      if (has_active_nodes && !has_mobile_nodes)
        {
          const unsigned int cell_id = (*cell)->active_cell_index();
          cell_mobility_status_map.insert({cell_id, mobility_status::active});
        }

      // Iterate to the next cell, none of the cells needs to be removed from
      // the cell set as done in previous loops since there are no other loop
      // over this set
      ++cell;
    }
}

template class DisableParticleContact<2>;
template class DisableParticleContact<3>;