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

template <int dim>
DisableParticleContact<dim>::DisableParticleContact()
{}

template <int dim>
void
DisableParticleContact<dim>::calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler)
{
  granular_temperature_average.clear();
  void_fractions.clear();
  granular_temperature_average.resize(triangulation.n_active_cells());
  void_fractions.resize(triangulation.n_active_cells());

  // Iterating through the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          double granular_temperature_cell = 0.0;
          double void_fraction             = 0.0;

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

              // Initialize variables for void fraction
              double       solid_volume = 0.0;
              const double cell_volume  = cell->measure();

              // Initialize velocity fluctuations
              Tensor<1, dim> cell_velocity_fluctuation_squared_sum;
              Tensor<1, dim> cell_velocity_fluctuation_squared_average;

              // First loop over particles in cell to calculation the average
              // velocity and the void fraction
              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_cell_iterator = particles_in_cell.begin();
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
                    M_PI *
                    pow(particle_properties[DEM::PropertiesIndex::dp], dim) /
                    (2.0 * dim);
                }

              // Calculate average velocity in the cell
              for (int d = 0; d < dim; ++d)
                velocity_cell_average[d] =
                  velocity_cell_sum[d] / n_particles_in_cell;

              // Calculate void fraction of cell
              void_fraction = 1 - solid_volume / cell_volume;

              // Second loop over particle to calculate the average granular
              // temperature
              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_cell_iterator = particles_in_cell.begin();
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
                    cell_velocity_fluctuation_squared_sum[d] /
                    n_particles_in_cell;
                  granular_temperature_cell +=
                    (1.0 / dim) * cell_velocity_fluctuation_squared_average[d];
                }
            }

          granular_temperature_average[cell->active_cell_index()] =
            granular_temperature_cell;
          void_fractions[cell->active_cell_index()] = void_fraction;
        }
    }
}

template <int dim>
void
DisableParticleContact<dim>::check_if_mobile(
  const typename dem_data_structures<dim>::cells_neighbor_list
    &                                    cells_neighbor_list,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto cell_neighbor_list_iterator = cells_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // 1st loop : assign mobility status to the main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      const unsigned int cell_id =
        (*cell_neighbor_iterator)->active_cell_index();
      const unsigned int n_particles_in_cell =
        particle_handler.n_particles_in_cell(*cell_neighbor_iterator);

      // Check to see if the cell has any particles
      if (n_particles_in_cell > 0)
        {
          if (granular_temperature_average[cell_id] >=
                granular_temperature_limit ||
              void_fractions[cell_id] >= void_fraction_limit)
            {
              cell_status[cell_id] = mobility_status::mobile;
              status_to_cell[mobility_status::mobile].insert(
                *cell_neighbor_iterator);
            }
        }

      if (cell_status[cell_id] != mobility_status::mobile)
        {
          status_to_cell[mobility_status::inactive].insert(
            *cell_neighbor_iterator);
        }

      // Assign mobility status to neighbor cells
      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          const unsigned int cell_id =
            (*cell_neighbor_iterator)->active_cell_index();

          const unsigned int n_particles_in_cell =
            particle_handler.n_particles_in_cell(*cell_neighbor_iterator);

          // Check if the cell has any particles
          if (n_particles_in_cell > 0)
            {
              if (granular_temperature_average[cell_id] >=
                    granular_temperature_limit ||
                  void_fractions[cell_id] >= void_fraction_limit)
                {
                  cell_status[cell_id] = mobility_status::mobile;
                  status_to_cell[mobility_status::mobile].insert(
                    *cell_neighbor_iterator);
                }
            }

          // No particles or not mobile => inactive
          if (cell_status[cell_id] != mobility_status::mobile)
            {
              status_to_cell[mobility_status::inactive].insert(
                *cell_neighbor_iterator);
            }
        }
    }
}


template <int dim>
void
DisableParticleContact<dim>::check_if_mobile_extended(
  const typename dem_data_structures<dim>::cells_neighbor_list
    &                                    cells_neighbor_list,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto cell_neighbor_list_iterator = cells_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      auto main_cell_iterator     = cell_neighbor_list_iterator->begin();
      unsigned int main_cell_id =
        (*cell_neighbor_iterator)->active_cell_index();
      const unsigned int n_particles_in_cell =
        particle_handler.n_particles_in_cell(*cell_neighbor_iterator);

      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          unsigned int neighbor_cell_id =
            (*cell_neighbor_iterator)->active_cell_index();

          const unsigned int n_particles_in_neighbor_cell =
            particle_handler.n_particles_in_cell(*cell_neighbor_iterator);

          // Cell is mobile if has particles and an empty neighbor
          if (n_particles_in_cell > 0 && n_particles_in_neighbor_cell == 0 &&
              cell_status[main_cell_id] != mobility_status::mobile)
            {
              cell_status[main_cell_id] = mobility_status::mobile;

              status_to_cell[mobility_status::inactive].erase(
                *main_cell_iterator);
              status_to_cell[mobility_status::mobile].insert(
                *main_cell_iterator);
            }


          //          if (cell_status[neighbor_cell_id] ==
          //          mobility_status::mobile &&
          //              cell_status[main_cell_id] ==
          //              mobility_status::inactive)
          //            {
          //              cell_status[main_cell_id] =
          //              mobility_status::mobile_extended;
          //
          //              status_to_cell[mobility_status::inactive].erase(
          //                *main_cell_iterator);
          //              status_to_cell[mobility_status::mobile_extended].insert(
          //                *main_cell_iterator);
          //            }

          if (n_particles_in_neighbor_cell > 0 && n_particles_in_cell == 0 &&
              cell_status[neighbor_cell_id] != mobility_status::mobile)
            {
              cell_status[neighbor_cell_id] = mobility_status::mobile;

              status_to_cell[mobility_status::inactive].erase(
                *cell_neighbor_iterator);
              status_to_cell[mobility_status::mobile].insert(
                *cell_neighbor_iterator);
            }

          //          // Neighbor is mobile if has particles and an empty cell
          //          if (cell_status[main_cell_id] == mobility_status::mobile
          //          &&
          //              cell_status[neighbor_cell_id] ==
          //              mobility_status::inactive)
          //            {
          //              cell_status[neighbor_cell_id] =
          //              mobility_status::mobile_extended;
          //
          //              status_to_cell[mobility_status::inactive].erase(
          //                *cell_neighbor_iterator);
          //              status_to_cell[mobility_status::mobile_extended].insert(
          //                *cell_neighbor_iterator);
          //            }
        }
    }
}

template <int dim>
void
DisableParticleContact<dim>::check_if_mobile_extended_extended(
  const typename dem_data_structures<dim>::cells_neighbor_list
    &                                    cells_neighbor_list,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  for (auto cell_neighbor_list_iterator = cells_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      auto main_cell_iterator     = cell_neighbor_list_iterator->begin();
      unsigned int main_cell_id =
        (*cell_neighbor_iterator)->active_cell_index();

      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          unsigned int neighbor_cell_id =
            (*cell_neighbor_iterator)->active_cell_index();

          if ((granular_temperature_average[main_cell_id] -
               granular_temperature_average[neighbor_cell_id]) >=
                0.25 * granular_temperature_limit &&
              cell_status[neighbor_cell_id] == mobility_status::inactive)
            {
              cell_status[neighbor_cell_id] = mobility_status::mobile;

              status_to_cell[mobility_status::inactive].erase(
                *cell_neighbor_iterator);
              status_to_cell[mobility_status::mobile].insert(
                *cell_neighbor_iterator);
            }

          if ((granular_temperature_average[neighbor_cell_id] -
               granular_temperature_average[main_cell_id]) >=
                0.25 * granular_temperature_limit &&
              cell_status[main_cell_id] == mobility_status::inactive)
            {
              cell_status[main_cell_id] = mobility_status::mobile;

              status_to_cell[mobility_status::inactive].erase(
                *main_cell_iterator);
              status_to_cell[mobility_status::mobile].insert(
                *main_cell_iterator);
            }
        }
    }
}

template <int dim>
void
DisableParticleContact<dim>::check_if_active(
  const typename dem_data_structures<dim>::cells_neighbor_list
    &cells_neighbor_list)
{
  for (auto cell_neighbor_list_iterator = cells_neighbor_list.begin();
       cell_neighbor_list_iterator != cells_neighbor_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      auto cell_neighbor_iterator = cell_neighbor_list_iterator->begin();
      auto main_cell_iterator     = cell_neighbor_list_iterator->begin();
      unsigned int main_cell_id =
        (*cell_neighbor_iterator)->active_cell_index();

      ++cell_neighbor_iterator;
      for (; cell_neighbor_iterator != cell_neighbor_list_iterator->end();
           ++cell_neighbor_iterator)
        {
          unsigned int neighbor_cell_id =
            (*cell_neighbor_iterator)->active_cell_index();

          if (cell_status[neighbor_cell_id] == mobility_status::mobile &&
              cell_status[main_cell_id] == mobility_status::inactive)
            {
              cell_status[main_cell_id] = mobility_status::active;

              status_to_cell[mobility_status::inactive].erase(
                *main_cell_iterator);
              status_to_cell[mobility_status::active].insert(
                *main_cell_iterator);
            }

          if (cell_status[main_cell_id] == mobility_status::mobile &&
              cell_status[neighbor_cell_id] == mobility_status::inactive)
            {
              cell_status[neighbor_cell_id] = mobility_status::active;

              status_to_cell[mobility_status::inactive].erase(
                *cell_neighbor_iterator);
              status_to_cell[mobility_status::active].insert(
                *cell_neighbor_iterator);
            }
        }
    }
}

template <int dim>
void
DisableParticleContact<dim>::identify_mobility_status(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler,
  const DEMContainerManager<dim> &                 container_manager)
{
  // Reset cell status containers
  cell_status.clear();
  cell_status.resize(triangulation.n_active_cells(), 0);

  status_to_cell.clear();
  status_to_cell.resize(mobility_status::n_mobility_status);

  std::list<typename dem_data_structures<dim>::cells_neighbor_list>
    all_neighbor_lists = {
      container_manager.cells_local_neighbor_list,
      container_manager.cells_ghost_neighbor_list,
      container_manager.cells_local_periodic_neighbor_list,
      container_manager.cells_ghost_periodic_neighbor_list,
      container_manager.cells_ghost_local_periodic_neighbor_list};

  calculate_average_granular_temperature(triangulation, particle_handler);

  // Check every status, one container at the time since each status check
  // depends on the last one
  for (auto &neighbor_list : all_neighbor_lists)
    check_if_mobile(neighbor_list, particle_handler);

  for (auto &neighbor_list : all_neighbor_lists)
    check_if_mobile_extended(neighbor_list, particle_handler);

  for (auto &neighbor_list : all_neighbor_lists)
    check_if_mobile_extended_extended(neighbor_list, particle_handler);

  for (auto &neighbor_list : all_neighbor_lists)
    check_if_active(neighbor_list);
}

template class DisableParticleContact<2>;
template class DisableParticleContact<3>;