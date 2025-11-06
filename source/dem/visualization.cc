// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/visualization.h>

#include <deal.II/numerics/data_out.h>

using namespace dealii;

template <int dim, typename PropertiesIndex>
Visualization<dim, PropertiesIndex>::Visualization() = default;

template <int dim, typename PropertiesIndex>
void
Visualization<dim, PropertiesIndex>::build_patches(
  dealii::Particles::ParticleHandler<dim>        &particle_handler,
  const std::vector<std::pair<std::string, int>> &properties)
{
  // We reserve 1 more name than the number of properties, since we will also be
  // writing the ID of the particles
  this->dataset_names.reserve(PropertiesIndex::n_properties + 1);

  // Adding ID to properties vector for visualization
  dataset_names.emplace_back("ID");

  // Defining property field position by iterating over properties.
  for (int field_position = 0; field_position < PropertiesIndex::n_properties;
       ++field_position)
    {
      // Get the property field name
      const std::string field_name = properties[field_position].first;

      // Number of components of the corresponding property
      const unsigned n_components = properties[field_position].second;

      // Check to see if the property is a vector
      // By default we assume that even 2D simulations have 3 components
      // Since the velocity and the angular velocity are stored in 3D
      // even for 2D simulations
      if (n_components == 3)
        {
          // The property is a vector, thus we set that the components
          // are part of a vector. Do not forget that since we added the ID
          // of the particles, we need to shift the field_position by 1. Hence,
          // we add 1 to field_position
          vector_datasets.emplace_back(
            field_position + 1,
            field_position + n_components,
            field_name,
            DataComponentInterpretation::component_is_part_of_vector);
        }
      dataset_names.emplace_back(field_name);
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
      patches[i].data.reinit(PropertiesIndex::n_properties + 1, 1);

      // ID and other properties
      if (particle->has_properties())
        {
          auto particle_properties = particle->get_properties();

          // Adding ID to patches
          patches[i].data(0, 0) = particle->get_id();

          // We need to offset the data we write by one since we have one more
          // property due to the fact that we save the ID of the particles.
          for (unsigned int property_index = 0;
               property_index < PropertiesIndex::n_properties;
               ++property_index)
            patches[i].data(property_index + 1, 0) =
              particle_properties[property_index];
        }
    }
}

template <int dim, typename PropertiesIndex>
void
Visualization<dim, PropertiesIndex>::print_xyz(
  dealii::Particles::ParticleHandler<dim> &particle_handler,
  const MPI_Comm                          &mpi_communicator,
  const ConditionalOStream                &pcout)
{
  const bool is_dem_mp =
    std::is_same_v<PropertiesIndex, DEM::DEMMPProperties::PropertiesIndex>;
  pcout << "id, type, dp, x, y, z";
  if constexpr (is_dem_mp)
    {
      pcout << ", T";
    }
  pcout << " " << std::endl;
  // Aggressively force synchronization of the header line
  usleep(500);
  MPI_Barrier(mpi_communicator);
  usleep(500);
  MPI_Barrier(mpi_communicator);

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
                        << particle_properties[PropertiesIndex::type] << " "
                        << std::setprecision(5)
                        << particle_properties[PropertiesIndex::dp] << " "
                        << std::setprecision(4) << particle_location << " ";
              if constexpr (is_dem_mp)
                {
                  std::cout << std::fixed << std::setprecision(4)
                            << particle_properties[PropertiesIndex::T];
                }
              std::cout << std::endl;
            }
        }
      usleep(500);
      MPI_Barrier(mpi_communicator);
    }
}

template <int dim, typename PropertiesIndex>
void
Visualization<dim, PropertiesIndex>::print_intermediate_format(
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

  // Add data in deal.II intermediate format to the string stream buffer for
  // each processor in order
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

template <int dim, typename PropertiesIndex>
const std::vector<DataOutBase::Patch<0, dim>> &
Visualization<dim, PropertiesIndex>::get_patches() const
{
  return patches;
}

template <int dim, typename PropertiesIndex>
std::vector<std::string>
Visualization<dim, PropertiesIndex>::get_dataset_names() const
{
  return dataset_names;
}

template <int dim, typename PropertiesIndex>
std::vector<
  std::tuple<unsigned int,
             unsigned int,
             std::string,
             DataComponentInterpretation::DataComponentInterpretation>>
Visualization<dim, PropertiesIndex>::get_nonscalar_data_ranges() const
{
  return vector_datasets;
}

template <int dim, typename PropertiesIndex>
Visualization<dim, PropertiesIndex>::~Visualization() = default;

template class Visualization<2, DEM::DEMProperties::PropertiesIndex>;
template class Visualization<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class Visualization<2, DEM::DEMMPProperties::PropertiesIndex>;
template class Visualization<3, DEM::DEMProperties::PropertiesIndex>;
template class Visualization<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class Visualization<3, DEM::DEMMPProperties::PropertiesIndex>;
