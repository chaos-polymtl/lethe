#include <core/solutions_output.h>

#include <dem/lagrangian_post_processing.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/block_vector.h>

#include <deal.II/numerics/data_out.h>

#include <fstream>

using namespace dealii;

template <int dim>
LagrangianPostProcessing<dim>::LagrangianPostProcessing()
{}

template <int dim>
void
LagrangianPostProcessing<dim>::calculate_average_particles_velocity(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler)
{
  velocity_average.reinit(triangulation.n_active_cells() * dim);

  // Iterating through the active cells in the triangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          Tensor<1, dim> velocity_in_cell_average =
            calculate_cell_average_particles_velocity(cell, particle_handler);
          for (unsigned int d = 0; d < dim; d++)
            {
              velocity_average[cell->active_cell_index() + d] =
                velocity_in_cell_average[d];
            }
        }
    }
}

template <int dim>
void
LagrangianPostProcessing<dim>::calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler)
{
  granular_temperature_average.reinit(triangulation.n_active_cells());

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
                calculate_cell_average_particles_velocity(cell,
                                                          particle_handler);

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
                  auto &particle_properties =
                    particles_in_cell_iterator->get_properties();

                  for (int d = 0; d < dim; ++d)
                    {
                      cell_velocity_fluctuation_squared_sum[d] +=
                        (particle_properties[DEM::PropertiesIndex::v_x + d] -
                         velocity_in_cell_average[d]) *
                        (particle_properties[DEM::PropertiesIndex::v_x + d] -
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

template <int dim>
Tensor<1, dim>
LagrangianPostProcessing<dim>::calculate_cell_average_particles_velocity(
  const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
  const Particles::ParticleHandler<dim> &particle_handler)
{
  Tensor<1, dim> velocity_cell_sum;
  Tensor<1, dim> velocity_cell_average;
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
          auto &particle_properties =
            particles_in_cell_iterator->get_properties();

          for (int d = 0; d < dim; ++d)
            {
              velocity_cell_sum[d] +=
                particle_properties[DEM::PropertiesIndex::v_x + d];
            }

          particles_cell_number++;
        }
      // Calculate average velocity in the cell
      for (int d = 0; d < dim; ++d)
        velocity_cell_average[d] = velocity_cell_sum[d] / particles_cell_number;
    }
  return velocity_cell_average;
}

template <int dim>
void
LagrangianPostProcessing<dim>::write_post_processing_results(
  const parallel::distributed::Triangulation<dim> &triangulation,
  PVDHandler &                                     grid_pvdhandler,
  const Particles::ParticleHandler<dim> &          particle_handler,
  const DEMSolverParameters<dim> &                 dem_parameters,
  DoFHandler<dim> &                                background_dh,
  const double                                     current_time,
  const unsigned int                               step_number,
  const MPI_Comm &                                 mpi_communicator)
{
  const std::string folder = dem_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    dem_parameters.simulation_control.output_name;
  const unsigned int group_files =
    dem_parameters.simulation_control.group_files;

  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);

  std::vector<std::string> average_solution_names(dim, "average_velocity");

  // Write particles' average velocity
  calculate_average_particles_velocity(triangulation, particle_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  data_out.add_data_vector(velocity_average,
                           average_solution_names,
                           DataOut<dim>::type_cell_data,
                           data_component_interpretation);


  // Write particles' granular temperature
  calculate_average_granular_temperature(triangulation, particle_handler);
  average_solution_names.push_back("granular_temperature");

  data_out.add_data_vector(granular_temperature_average,
                           average_solution_names.back(),
                           DataOut<dim>::type_cell_data);

  // Write grid
  data_out.attach_dof_handler(background_dh);

  // Attach the solution data to data_out object
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");


  const std::string postprocess_file_name =
    dem_parameters.simulation_control.output_name + "-postprocess_data";

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

template class LagrangianPostProcessing<2>;
template class LagrangianPostProcessing<3>;
