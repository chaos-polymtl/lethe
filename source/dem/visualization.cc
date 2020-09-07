/*
 * Visualization.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include <dem/visualization.h>

using namespace dealii;

template <int dim>
Visualization<dim>::Visualization()
{}

template <int dim>
void
Visualization<dim>::build_patches(
  const dealii::Particles::ParticleHandler<dim> &particle_handler,
  std::vector<std::pair<std::string, int>>       properties)
{
  // Get the number of properties as a local parameter
  const int properties_number = particle_handler.n_properties_per_particle();

  // Defining property field position
  int field_position = 0;

  // Iterating over properties
  for (auto properties_iterator = properties.begin();
       properties_iterator != properties.end();
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
      patches[i].vertices[0]    = particle->get_location();
      patches[i].patch_index    = i;
      patches[i].n_subdivisions = 1;
      patches[i].data.reinit(properties_number, 1);

      // Other properties
      if (particle->has_properties())
        {
          const ArrayView<const double> props = particle->get_properties();

          for (unsigned int property_index = 0; property_index < props.size();
               ++property_index)
            patches[i].data(property_index, 0) = props[property_index];
        }
    }
}

template <int dim>
void
Visualization<dim>::print_xyz(
  const dealii::Particles::ParticleHandler<dim> &particle_handler,
  std::vector<std::pair<std::string, int>>       properties)
{
  std::vector<int> precision = {
    0, 0, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3,
  };
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    std::cout
      << "id , type, dp   , rho  , v_x  , v_y  , v_z  , acc_x , acc_y , "
         "acc_z , force_x, force_y, force_z, omega_x, omega_y, omega_z, "
         " mass , mom_inertia, M_x  , M_y  , M_z  "
      << std::endl;

  this->build_patches(particle_handler, properties);
  unsigned int counter;

  // loop over all patches
  for (const auto &patch : patches)
    {
      counter                  = 0;
      unsigned int n_data_sets = properties.size();
      for (unsigned int data_set = 0; data_set < n_data_sets;
           ++data_set, ++counter)
        {
          std::cout.precision(precision[counter]);
          std::cout << std::fixed << patch.data(data_set, 0) << " ";
        }
      std::cout << '\n';
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
