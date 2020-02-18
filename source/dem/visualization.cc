/*
 * Visualization.cpp
 *
 *  Created on: Oct 1, 2019
 *      Author: shahab
 */

#include "dem/visualization.h"

using namespace dealii;

template <int dim, int spacedim>
Visualization<dim, spacedim>::Visualization() {}

template <int dim, int spacedim>
void Visualization<dim, spacedim>::build_patches(
    const dealii::Particles::ParticleHandler<dim, spacedim> &particle_handler,
    std::vector<std::pair<std::string, int>> properties) {

  // Get the number of properties as a local parameter
  const int properties_number = particle_handler.n_properties_per_particle();

  // Defining property field position
  int field_position = 0;

  // Iterating over properties
  for (auto properties_iterator = properties.begin();
       properties_iterator != properties.end();
       ++properties_iterator, ++field_position) {

    // Get the property field name
    const std::string field_name = properties_iterator->first;

    // Number of components of the corresponding property
    const unsigned components_number = properties_iterator->second;

    // Check to see if the property is a vector
    if (components_number == dim) {
      vector_datasets.push_back(std::make_tuple(
          field_position, field_position + components_number - 1, field_name,
          DataComponentInterpretation::component_is_part_of_vector));
    }
    dataset_names.push_back(field_name);
  }

  // Building the patch data
  patches.resize(particle_handler.n_locally_owned_particles());

  typename dealii::Particles::ParticleHandler<dim, spacedim>::particle_iterator
      particle = particle_handler.begin();

  // Looping over particle to get the properties from the particle_handler
  for (unsigned int i = 0; particle != particle_handler.end();
       ++particle, ++i) {
    // Particle location
    patches[i].vertices[0] = particle->get_location();
    patches[i].patch_index = i;
    patches[i].n_subdivisions = 1;
    patches[i].data.reinit(properties_number, 1);

    // Other properties
    if (particle->has_properties()) {
      const ArrayView<const double> props = particle->get_properties();

      for (unsigned int property_index = 0; property_index < props.size();
           ++property_index)
        patches[i].data(property_index, 0) = props[property_index];
    }
  }
}

template <int dim, int spacedim>
const std::vector<DataOutBase::Patch<0, dim>> &
Visualization<dim, spacedim>::get_patches() const {
  return patches;
}

template <int dim, int spacedim>
std::vector<std::string>
Visualization<dim, spacedim>::get_dataset_names() const {
  return dataset_names;
}

template <int dim, int spacedim>
std::vector<
    std::tuple<unsigned int, unsigned int, std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
Visualization<dim, spacedim>::get_nonscalar_data_ranges() const {
  return vector_datasets;
}

template <int dim, int spacedim>
Visualization<dim, spacedim>::~Visualization() {}

template class Visualization<3, 3>;
