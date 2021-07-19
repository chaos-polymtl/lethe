#include <dem/lagrangian_post_processing.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/numerics/data_out.h>

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
  Vector<double> velocity_average(dim * triangulation.n_active_cells());

  // Iterating through the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          Tensor<1, dim> cell_velocity_average =
            calculate_cell_average_particles_velocity(cell, particle_handler);

          for (int d = 0; d < dim; ++d)
            velocity_average[dim * cell->global_active_cell_index() + d] =
              cell_velocity_average[d];
        }
    }

  // Writing data
  DataOut<dim> data_out;
  data_out.attach_triangulation(triangulation);

  std::vector<std::string> average_solution_names(dim, "average_velocity");

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    average_data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  data_out.add_data_vector(velocity_average,
                           average_solution_names,
                           DataOut<dim>::type_cell_data,
                           average_data_component_interpretation);


  // data_out.build_patches();
  // std::ofstream output("solution-" + std::to_string(10) + ".vtk");
  // data_out.write_vtk(output);
}

template <int dim>
void
LagrangianPostProcessing<dim>::calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler)
{
  Vector<double> granular_temperature_average(triangulation.n_cells());

  double       granular_temperature_cell(0);
  unsigned int particles_cell_number(0);

  // Iterating through the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
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
                    (1 / dim) * cell_velocity_fluctuation_squared_average[d];
                }
            }
          granular_temperature_average[cell->global_active_cell_index()] =
            granular_temperature_cell;
        }
    }

  // Writing output file
}

template <int dim>
Tensor<1, dim>
LagrangianPostProcessing<dim>::calculate_cell_average_particles_velocity(
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

template class LagrangianPostProcessing<2>;
template class LagrangianPostProcessing<3>;
