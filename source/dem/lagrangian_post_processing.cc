// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>
#include <core/solutions_output.h>

#include <dem/lagrangian_post_processing.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
void
calculate_average_particles_velocity(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim>           &particle_handler,
  Vector<double>                                  &velocity_average_x,
  Vector<double>                                  &velocity_average_y,
  Vector<double>                                  &velocity_average_z,
  Vector<double>                                  &velocity_average_magnitude)
{
  // Iterating through the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          Tensor<1, dim> cell_velocity_average =
            calculate_cell_average_particles_velocity<dim, PropertiesIndex>(
              cell, particle_handler);

          velocity_average_x[cell->active_cell_index()] =
            cell_velocity_average[0];
          velocity_average_y[cell->active_cell_index()] =
            cell_velocity_average[1];

          if constexpr (dim == 3)
            {
              velocity_average_z[cell->active_cell_index()] =
                cell_velocity_average[2];
              velocity_average_magnitude[cell->active_cell_index()] =
                sqrt(pow(cell_velocity_average[0], 2) +
                     pow(cell_velocity_average[1], 2) +
                     pow(cell_velocity_average[2], 2));
            }

          if constexpr (dim == 2)
            velocity_average_magnitude[cell->active_cell_index()] =
              sqrt(pow(cell_velocity_average[0], 2) +
                   pow(cell_velocity_average[1], 2));
        }
    }
}

template <int dim, typename PropertiesIndex>
void
calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim>           &particle_handler,
  Vector<double>                                  &granular_temperature_average)
{
  // Iterating through the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          double       granular_temperature_cell(0);
          unsigned int particles_cell_number(0);

          // Looping through all the particles in the cell
          // Particles in the cell
          typename Particles::ParticleHandler<dim>::particle_iterator_range
            particles_in_cell = particle_handler.particles_in_cell(cell);

          const bool particles_exist_in_main_cell = !particles_in_cell.empty();
          // Check to see if the cell has any particles
          if (particles_exist_in_main_cell)
            {
              Tensor<1, dim> velocity_in_cell_average =
                calculate_cell_average_particles_velocity<dim, PropertiesIndex>(
                  cell, particle_handler);

              // Initializing velocity fluctuations
              Tensor<1, dim> cell_velocity_fluctuation_squared_sum =
                Tensor<1, dim>();
              Tensor<1, dim> cell_velocity_fluctuation_squared_average =
                Tensor<1, dim>();

              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_cell_iterator = particles_in_cell.begin();
                   particles_in_cell_iterator != particles_in_cell.end();
                   ++particles_in_cell_iterator)
                {
                  auto particle_properties =
                    particles_in_cell_iterator->get_properties();

                  for (int d = 0; d < dim; ++d)
                    {
                      cell_velocity_fluctuation_squared_sum[d] +=
                        (particle_properties[PropertiesIndex::v_x + d] -
                         velocity_in_cell_average[d]) *
                        (particle_properties[PropertiesIndex::v_x + d] -
                         velocity_in_cell_average[d]);
                    }

                  particles_cell_number++;
                }
              // Calculate average granular temperature in the cell
              for (int d = 0; d < dim; ++d)
                {
                  cell_velocity_fluctuation_squared_average[d] =
                    cell_velocity_fluctuation_squared_sum[d] /
                    particles_cell_number;
                  granular_temperature_cell +=
                    (1.0 / dim) * cell_velocity_fluctuation_squared_average[d];
                }
            }
          granular_temperature_average[cell->active_cell_index()] =
            granular_temperature_cell;
        }
    }
}

template <int dim, typename PropertiesIndex>
Tensor<1, dim>
calculate_cell_average_particles_velocity(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  Tensor<1, dim> velocity_cell_sum     = Tensor<1, dim>();
  Tensor<1, dim> velocity_cell_average = Tensor<1, dim>();
  unsigned int   particles_cell_number(0);

  // Looping through all the particles in the cell
  // Particles in the cell
  typename Particles::ParticleHandler<dim>::particle_iterator_range
    particles_in_cell = particle_handler.particles_in_cell(cell);

  const bool particles_exist_in_main_cell = !particles_in_cell.empty();

  // Check to see if the cell has any particles
  if (particles_exist_in_main_cell)
    {
      for (typename Particles::ParticleHandler<dim>::particle_iterator_range::
             iterator particles_in_cell_iterator = particles_in_cell.begin();
           particles_in_cell_iterator != particles_in_cell.end();
           ++particles_in_cell_iterator)
        {
          auto particle_properties =
            particles_in_cell_iterator->get_properties();

          for (int d = 0; d < dim; ++d)
            {
              velocity_cell_sum[d] +=
                particle_properties[PropertiesIndex::v_x + d];
            }

          particles_cell_number++;
        }
      // Calculate average velocity in the cell
      for (int d = 0; d < dim; ++d)
        velocity_cell_average[d] = velocity_cell_sum[d] / particles_cell_number;
    }
  return velocity_cell_average;
}

template <int dim, typename PropertiesIndex>
void
write_post_processing_results(
  const parallel::distributed::Triangulation<dim> &triangulation,
  PVDHandler                                      &grid_pvdhandler,
  const DoFHandler<dim>                           &background_dh,
  const Particles::ParticleHandler<dim>           &particle_handler,
  const DEMSolverParameters<dim>                  &dem_parameters,
  const double                                     current_time,
  const unsigned int                               step_number,
  const MPI_Comm                                  &mpi_communicator,
  AdaptiveSparseContacts<dim, PropertiesIndex>    &sparse_contacts_object)
{
  const std::string folder = dem_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    dem_parameters.simulation_control.output_name;
  const unsigned int group_files =
    dem_parameters.simulation_control.group_files;

  DataOut<dim> data_out;
  data_out.attach_dof_handler(background_dh);

  // Write particles' average velocity
  Vector<double> velocity_average_x(triangulation.n_active_cells());
  Vector<double> velocity_average_y(triangulation.n_active_cells());
  Vector<double> velocity_average_z(triangulation.n_active_cells());
  Vector<double> velocity_average_magnitude(triangulation.n_active_cells());
  calculate_average_particles_velocity<dim, PropertiesIndex>(
    triangulation,
    particle_handler,
    velocity_average_x,
    velocity_average_y,
    velocity_average_z,
    velocity_average_magnitude);

  data_out.add_data_vector(velocity_average_x,
                           "average_velocity_x",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(velocity_average_y,
                           "average_velocity_y",
                           DataOut<dim>::type_cell_data);

  if constexpr (dim == 3)
    data_out.add_data_vector(velocity_average_z,
                             "average_velocity_z",
                             DataOut<dim>::type_cell_data);

  data_out.add_data_vector(velocity_average_magnitude,
                           "average_velocity_magnitude",
                           DataOut<dim>::type_cell_data);

  // Write particles' granular temperature
  Vector<double> granular_temperature_average(triangulation.n_active_cells());
  calculate_average_granular_temperature<dim, PropertiesIndex>(
    triangulation, particle_handler, granular_temperature_average);


  data_out.add_data_vector(granular_temperature_average,
                           "granular_temperature",
                           DataOut<dim>::type_cell_data);

  // Write mobility status of cells
  Vector<float> mobility_status(triangulation.n_active_cells());
  sparse_contacts_object.get_mobility_status_vector(mobility_status);
  data_out.add_data_vector(mobility_status,
                           "mobility_status",
                           DataOut<dim>::type_cell_data);

  // Attach the solution data to data_out object
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");

  const std::string postprocess_file_name =
    dem_parameters.simulation_control.output_name +
    "_lagrangian_postprocessing";

  data_out.build_patches();

  write_vtu_and_pvd<dim>(grid_pvdhandler,
                         data_out,
                         folder,
                         postprocess_file_name,
                         current_time,
                         step_number,
                         group_files,
                         mpi_communicator);
}

template void
write_post_processing_results<2, DEM::DEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<2>                           &background_dh,
  const Particles::ParticleHandler<2>           &particle_handler,
  const DEMSolverParameters<2>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<2, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
write_post_processing_results<2, DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<2>                           &background_dh,
  const Particles::ParticleHandler<2>           &particle_handler,
  const DEMSolverParameters<2>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<2, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
write_post_processing_results<2, DEM::DEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<2>                           &background_dh,
  const Particles::ParticleHandler<2>           &particle_handler,
  const DEMSolverParameters<2>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<2, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
write_post_processing_results<3, DEM::DEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<3>                           &background_dh,
  const Particles::ParticleHandler<3>           &particle_handler,
  const DEMSolverParameters<3>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<3, DEM::DEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
write_post_processing_results<3, DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<3>                           &background_dh,
  const Particles::ParticleHandler<3>           &particle_handler,
  const DEMSolverParameters<3>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<3, DEM::CFDDEMProperties::PropertiesIndex>
    &sparse_contacts_object);

template void
write_post_processing_results<3, DEM::DEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3> &triangulation,
  PVDHandler                                    &grid_pvdhandler,
  const DoFHandler<3>                           &background_dh,
  const Particles::ParticleHandler<3>           &particle_handler,
  const DEMSolverParameters<3>                  &dem_parameters,
  const double                                   current_time,
  const unsigned int                             step_number,
  const MPI_Comm                                &mpi_communicator,
  AdaptiveSparseContacts<3, DEM::DEMMPProperties::PropertiesIndex>
    &sparse_contacts_object);

template Tensor<1, 2>
calculate_cell_average_particles_velocity<2,
                                          DEM::DEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2>::cell_iterator &cell,
  const Particles::ParticleHandler<2> &particle_handler);

template Tensor<1, 2>
calculate_cell_average_particles_velocity<
  2,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2>::cell_iterator &cell,
  const Particles::ParticleHandler<2> &particle_handler);

template Tensor<1, 2>
calculate_cell_average_particles_velocity<
  2,
  DEM::DEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<2>::cell_iterator &cell,
  const Particles::ParticleHandler<2> &particle_handler);

template Tensor<1, 3>
calculate_cell_average_particles_velocity<3,
                                          DEM::DEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3>::cell_iterator &cell,
  const Particles::ParticleHandler<3> &particle_handler);

template Tensor<1, 3>
calculate_cell_average_particles_velocity<
  3,
  DEM::CFDDEMProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3>::cell_iterator &cell,
  const Particles::ParticleHandler<3> &particle_handler);

template Tensor<1, 3>
calculate_cell_average_particles_velocity<
  3,
  DEM::DEMMPProperties::PropertiesIndex>(
  const parallel::distributed::Triangulation<3>::cell_iterator &cell,
  const Particles::ParticleHandler<3> &particle_handler);
