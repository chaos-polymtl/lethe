// SPDX-FileCopyrightText: Copyright (c) 2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/adaptive_sparse_contacts.h>
#include <dem/dem_action_manager.h>
#include <dem/load_balancing.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
LagrangianLoadBalancing<dim, PropertiesIndex>::LagrangianLoadBalancing()
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
{}

template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::check_load_balance_once() const
{
  if (simulation_control->get_step_number() == load_balance_step)
    DEMActionManager::get_action_manager()->load_balance_step();
}

template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::check_load_balance_frequent()
  const
{
  if ((simulation_control->get_step_number() % load_balance_frequency) == 0)
    DEMActionManager::get_action_manager()->load_balance_step();
}

template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::check_load_balance_dynamic()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency == 0)
    {
      unsigned int maximum_particle_number_on_proc =
        Utilities::MPI::max(particle_handler->n_locally_owned_particles(),
                            mpi_communicator);
      unsigned int minimum_particle_number_on_proc =
        Utilities::MPI::min(particle_handler->n_locally_owned_particles(),
                            mpi_communicator);
      unsigned int average_particle_number_on_proc =
        particle_handler->n_locally_owned_particles() / n_mpi_processes;

      // Execute load balancing if difference of load between processors is
      // larger than threshold of the load per processor
      if ((maximum_particle_number_on_proc - minimum_particle_number_on_proc) >
          load_threshold * average_particle_number_on_proc)
        DEMActionManager::get_action_manager()->load_balance_step();
    }
}

template <int dim, typename PropertiesIndex>
inline void
LagrangianLoadBalancing<dim, PropertiesIndex>::
  check_load_balance_with_sparse_contacts()
{
  if (simulation_control->get_step_number() % dynamic_check_frequency == 0)
    {
      // Process to accumulate the load of each process regards the number
      // of cells and particles with their selected weight and with a factor
      // related to the mobility status of the cells
      double load_weight = 0.0;

      for (const auto &cell : triangulation->active_cell_iterators())
        {
          if (cell->is_locally_owned())
            {
              // Apply the cell weight
#if DEAL_II_VERSION_GTE(9, 7, 0)
              Point<dim> cell_barycenter = cell->center();
              load_weight +=
                static_cast<int>(cell_weight_function->value(cell_barycenter));
#else
              load_weight += cell_weight;
#endif


              // Get the mobility status of the cell and the number of particles
              const unsigned int cell_mobility_status =
                adaptive_sparse_contacts->check_cell_mobility(cell);
              const unsigned int n_particles_in_cell =
                particle_handler->n_particles_in_cell(cell);

              // Apply a factor on the particle weight regards the
              // mobility status. alpha = 1 by default for mobile cell, but
              // is modified if cell is active or inactive
              double alpha = 1.0;
              if (cell_mobility_status ==
                    AdaptiveSparseContacts<dim,
                                           PropertiesIndex>::static_active ||
                  cell_mobility_status ==
                    AdaptiveSparseContacts<dim,
                                           PropertiesIndex>::advected_active)
                {
                  alpha = active_status_factor;
                }
              else if (cell_mobility_status ==
                         AdaptiveSparseContacts<dim,
                                                PropertiesIndex>::inactive ||
                       cell_mobility_status ==
                         AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
                {
                  alpha = inactive_status_factor;
                }

              // Add the particle weight time the number of particles in the
              // cell to the processor load
              load_weight += alpha * n_particles_in_cell * particle_weight;
            }
        }

      // Find the minimum load on a processor
      double maximum_load_on_proc =
        Utilities::MPI::max(load_weight, mpi_communicator);

      // Find the minimum load on a processor
      double minimum_load_on_proc =
        Utilities::MPI::min(load_weight, mpi_communicator);

      // Get the total load
      double total_load = Utilities::MPI::sum(load_weight, mpi_communicator);

      if ((maximum_load_on_proc - minimum_load_on_proc) >
          load_threshold * (total_load / n_mpi_processes))
        {
          // Clear and connect a new cell weight function
          connect_mobility_status_weight_signals();

          DEMActionManager::get_action_manager()->load_balance_step();
        }
    }

  // Clear and connect a new cell weight function with new mobility status
  connect_mobility_status_weight_signals();
}

template <int dim, typename PropertiesIndex>
unsigned int
LagrangianLoadBalancing<dim, PropertiesIndex>::calculate_total_cell_weight(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const CellStatus status) const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  switch (status)
    {
      case CellStatus::cell_will_persist:
      case CellStatus::cell_will_be_refined:
        // If CELL_PERSIST, do as CELL_REFINE
        {
          const unsigned int n_particles_in_cell =
            particle_handler->n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
        }
      case CellStatus::cell_invalid:
        break;

      case CellStatus::children_will_be_coarsened:
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler->n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
        }
      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

template <int dim, typename PropertiesIndex>
unsigned int
LagrangianLoadBalancing<dim, PropertiesIndex>::
  calculate_total_cell_weight_with_mobility_status(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
                    &cell,
    const CellStatus status) const
{
  // Assign no weight to cells we do not own.
  if (!cell->is_locally_owned())
    return 0;

  // Get mobility status of the cell
  const unsigned int cell_mobility_status =
    adaptive_sparse_contacts->check_cell_mobility(cell);

  // Applied a factor on the particle weight regards the mobility status
  // Factor of 1 when mobile cell
  double alpha = 1.0;
  if (cell_mobility_status ==
        AdaptiveSparseContacts<dim, PropertiesIndex>::static_active ||
      cell_mobility_status ==
        AdaptiveSparseContacts<dim, PropertiesIndex>::advected_active)
    {
      alpha = active_status_factor;
    }
  else if (cell_mobility_status ==
             AdaptiveSparseContacts<dim, PropertiesIndex>::inactive ||
           cell_mobility_status ==
             AdaptiveSparseContacts<dim, PropertiesIndex>::advected)
    {
      alpha = inactive_status_factor;
    }

  switch (status)
    {
      case dealii::CellStatus::cell_will_persist:
      case dealii::CellStatus::cell_will_be_refined:
        {
          const unsigned int n_particles_in_cell =
            particle_handler->n_particles_in_cell(cell);
          return alpha * n_particles_in_cell * particle_weight;
        }
      case dealii::CellStatus::cell_invalid:
        break;
      case dealii::CellStatus::children_will_be_coarsened:
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0;
               child_index < GeometryInfo<dim>::max_children_per_cell;
               ++child_index)
            n_particles_in_cell +=
              particle_handler->n_particles_in_cell(cell->child(child_index));

          return alpha * n_particles_in_cell * particle_weight;
        }
      default:
        Assert(false, ExcInternalError());
        break;
    }

  return 0;
}

template class LagrangianLoadBalancing<2, DEM::DEMProperties::PropertiesIndex>;
template class LagrangianLoadBalancing<2,
                                       DEM::CFDDEMProperties::PropertiesIndex>;
template class LagrangianLoadBalancing<2,
                                       DEM::DEMMPProperties::PropertiesIndex>;
template class LagrangianLoadBalancing<3, DEM::DEMProperties::PropertiesIndex>;
template class LagrangianLoadBalancing<3,
                                       DEM::CFDDEMProperties::PropertiesIndex>;
template class LagrangianLoadBalancing<3,
                                       DEM::DEMMPProperties::PropertiesIndex>;
