#include <dem/visualization.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

template <int dim>
Visualization<dim>::Visualization()
{}

template <int dim>
void
Visualization<dim>::build_patches(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  std::vector<std::pair<std::string, int>> properties)
{
  // Adding ID to properties vector for visualization
  properties.insert(properties.begin(), std::make_pair("ID", 1));

  // Defining properties for writing
  this->properties_to_write.assign(properties.begin(),
                                   properties.begin() +
                                     DEM::get_number_properties());

  // Defining property field position
  int field_position = 0;
  // Iterating over properties
  for (auto properties_iterator = properties_to_write.begin();
       properties_iterator != properties_to_write.end();
       ++properties_iterator, ++field_position)
    {
      // Get the property field name
      const std::string field_name = properties_iterator->first;

      // Number of components of the corresponding property
      const unsigned components_number = properties_iterator->second;

      // Check to see if the property is a vector
      if (components_number == dim)
        {
          vector_datasets.push_back(std::make_tuple(
            field_position,
            field_position + components_number - 1,
            field_name,
            DataComponentInterpretation::component_is_part_of_vector));
        }
      dataset_names.push_back(field_name);
    }

  // Building the patch data
  patches.resize(particle_handler.n_locally_owned_particles());

  typename dealii::Particles::ParticleHandler<dim>::particle_iterator particle =
    particle_handler.begin();

  // Looping over particle to get the properties from the particle_handler
  for (unsigned int i = 0; particle != particle_handler.end(); ++particle, ++i)
    {
      // Particle location
      patches[i].vertices[0] = particle->get_location();
      patches[i].patch_index = i;
      patches[i].data.reinit(DEM::get_number_properties(), 1);

      // ID and other properties
      if (particle->has_properties())
        {
          // Calculating force for visualization
          auto particle_properties = particle->get_properties();

          // Adding ID to patches
          patches[i].data(0, 0) = particle->get_id();

          for (unsigned int property_index = 1;
               property_index < DEM::get_number_properties();
               ++property_index)
            patches[i].data(property_index, 0) =
              particle_properties[property_index - 1];
        }
    }
}

template <int dim>
void
Visualization<dim>::print_xyz(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const MPI_Comm                          &mpi_communicator,
  const ConditionalOStream                &pcout)
{
  pcout << "id, type, dp, x, y, z " << std::endl;
  usleep(100);

  std::map<int, Particles::ParticleIterator<dim>> global_particles;
  unsigned int current_id, current_id_max = 0;

  // Mapping of all particles & find the max id on current processor
  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      current_id     = particle->get_id();
      current_id_max = std::max(current_id, current_id_max);

      global_particles.insert({current_id, particle});
    }

  // Find global max particle index
  unsigned int id_max = Utilities::MPI::max(current_id_max, mpi_communicator);

  // Print particle info one by one in ascending order
  for (unsigned int i = 0; i <= id_max; i++)
    {
      for (auto &iterator : global_particles)
        {
          unsigned int id = iterator.first;
          if (id == i)
            {
              auto particle            = iterator.second;
              auto particle_properties = particle->get_properties();
              auto particle_location   = particle->get_location();

              std::cout << std::fixed << std::setprecision(0) << id << " "
                        << std::setprecision(0)
                        << particle_properties[DEM::PropertiesIndex::type]
                        << " " << std::setprecision(5)
                        << particle_properties[DEM::PropertiesIndex::dp] << " "
                        << std::setprecision(4) << particle_location
                        << std::endl;
            }
        }
      usleep(100);
      MPI_Barrier(mpi_communicator);
    }
}

template <int dim>
void
Visualization<dim>::print_intermediate_format(
  const Vector<float>   &data_to_print,
  const DoFHandler<dim> &background_dh,
  const MPI_Comm        &mpi_communicator)
{
  unsigned int n_mpi_processes(
    Utilities::MPI::n_mpi_processes(mpi_communicator));
  unsigned int this_mpi_process(
    Utilities::MPI::this_mpi_process(mpi_communicator));

  // Attach the data to the background mesh
  DataOut<dim> data_out;
  data_out.attach_dof_handler(background_dh);
  data_out.add_data_vector(data_to_print,
                           "print_from_processor_" +
                             Utilities::int_to_string(this_mpi_process),
                           DataOut<dim>::type_cell_data);
  data_out.build_patches();

  std::stringstream out;

  // Add data in deal.II intermediate format to the stringstream buffer for each
  // processor in order
  for (unsigned int processor_number = 0; processor_number < n_mpi_processes;
       ++processor_number)
    {
      usleep(100);
      MPI_Barrier(mpi_communicator);
      if (processor_number == this_mpi_process)
        {
          // Generate a string stream buffer to store the output
          data_out.write_deal_II_intermediate(out);

          // Print in terminal but remove part of the header since it
          // contains some deal.II version information that may change
          std::string  line;
          unsigned int counter = 0;
          while (std::getline(out, line))
            {
              if (counter++ > 4)
                std::cout << line << std::endl;
            }
        }
    }
}

template <int dim>
const std::vector<DataOutBase::Patch<0, dim>> &
Visualization<dim>::get_patches() const
{
  return patches;
}

template <int dim>
std::vector<std::string>
Visualization<dim>::get_dataset_names() const
{
  return dataset_names;
}

template <int dim>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
Visualization<dim>::get_nonscalar_data_ranges() const
{
  return vector_datasets;
}

template <int dim>
Visualization<dim>::~Visualization()
{}

template class Visualization<2>;
template class Visualization<3>;
