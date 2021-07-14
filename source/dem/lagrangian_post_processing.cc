#include <dem/dem_properties.h>
#include <dem/lagrangian_post_processing.h>


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
  std::vector<Tensor<1, dim>> velocity_average;
  velocity_average.reserve(triangulation.n_cells());

  // Iterating through the active cells in the trangulation
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        velocity_average[cell->global_active_cell_index()] =
          calculate_cell_average_particles_velocity(cell, particle_handler);
    }

  // Writing output file
}

template <int dim>
void
LagrangianPostProcessing<dim>::calculate_average_granular_temperature(
  const parallel::distributed::Triangulation<dim> &triangulation,
  const Particles::ParticleHandler<dim> &          particle_handler)
{
  std::vector<double> granular_temperature_average;
  granular_temperature_average.reserve(triangulation.n_cells());

  double       granular_temperature_cell = 0;
  unsigned int particles_cell_number     = 0;

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
              Tensor<1, dim> cell_velocity_fluctuation_squared_sum;
              cell_velocity_fluctuation_squared_sum = 0;
              Tensor<1, dim> cell_velocity_fluctuation_squared_average;
              cell_velocity_fluctuation_squared_average = 0;

              for (typename Particles::ParticleHandler<
                     dim>::particle_iterator_range::iterator
                     particles_in_cell_iterator = particles_in_cell.begin();
                   particles_in_cell_iterator != particles_in_cell.end();
                   ++particles_in_cell_iterator)
                {
                  auto &particle_properties =
                    particles_in_cell_iterator->get_properties();

                  cell_velocity_fluctuation_squared_sum[0] +=
                    (particle_properties[DEM::PropertiesIndex::v_x] -
                     velocity_in_cell_average[0]) *
                    (particle_properties[DEM::PropertiesIndex::v_x] -
                     velocity_in_cell_average[0]);
                  cell_velocity_fluctuation_squared_sum[1] +=
                    (particle_properties[DEM::PropertiesIndex::v_y] -
                     velocity_in_cell_average[1]) *
                    (particle_properties[DEM::PropertiesIndex::v_y] -
                     velocity_in_cell_average[1]);
                  if (dim == 3)
                    cell_velocity_fluctuation_squared_sum[2] +=
                      (particle_properties[DEM::PropertiesIndex::v_z] -
                       velocity_in_cell_average[2]) *
                      (particle_properties[DEM::PropertiesIndex::v_z] -
                       velocity_in_cell_average[2]);

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
  Tensor<1, dim> velocity_cell_sum;
  velocity_cell_sum = 0;
  Tensor<1, dim> velocity_cell_average;
  velocity_cell_average              = 0;
  unsigned int particles_cell_number = 0;

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

          velocity_cell_sum[0] +=
            particle_properties[DEM::PropertiesIndex::v_x];
          velocity_cell_sum[1] +=
            particle_properties[DEM::PropertiesIndex::v_y];
          if (dim == 3)
            velocity_cell_sum[2] +=
              particle_properties[DEM::PropertiesIndex::v_z];

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
