// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/insertion_file.h>

using namespace DEM;

template <int dim, typename PropertiesIndex>
InsertionFile<dim, PropertiesIndex>::InsertionFile(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
  : Insertion<dim, PropertiesIndex>(size_distribution_object_container,
                                    triangulation,
                                    dem_parameters)
  , remaining_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
  , number_of_files(dem_parameters.insertion_info.list_of_input_files.size())
  , insertion_files(dem_parameters.insertion_info.list_of_input_files)
{
  // Initializing current inserting particle type and file id
  this->current_inserting_particle_type = 0;
  this->current_file_id                 = 0;
}

template <int dim, typename PropertiesIndex>
void
InsertionFile<dim, PropertiesIndex>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  if (remaining_particles_of_each_type == 0 &&
      this->current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remaining_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++this->current_inserting_particle_type);
    }

  if (remaining_particles_of_each_type > 0)
    {
      if (this->removing_particles_in_region)
        {
          if (this->mark_for_update)
            {
              this->find_cells_in_removing_box(triangulation);
              this->mark_for_update = false;
            }
          this->remove_particles_in_box(particle_handler);
        }

      MPI_Comm  communicator     = triangulation.get_mpi_communicator();
      const int this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

      // Read the input file on Rank 0
      std::map<std::string, std::vector<double>> particles_data;
      if (this_mpi_process == 0)
        {
          fill_vectors_from_file(particles_data,
                                 insertion_files.at(current_file_id),
                                 ";");
        }

      // Increment file id for next call
      current_file_id++;
      current_file_id = current_file_id % number_of_files;

      // Broadcast keys to all processes to reconstruct the map structure
      unsigned int n_keys = 0;
      if (this_mpi_process == 0)
        {
          n_keys = particles_data.size();
        }
      MPI_Bcast(&n_keys, 1, MPI_UNSIGNED, 0, communicator);

      std::vector<std::string> keys;
      if (this_mpi_process == 0)
        {
          for (auto const &[key, val] : particles_data)
            {
              keys.push_back(key);
            }
        }
      else
        {
          keys.resize(n_keys);
        }

      for (unsigned int i = 0; i < n_keys; ++i)
        {
          unsigned int key_size = 0;
          if (this_mpi_process == 0)
            {
              key_size = keys[i].size();
            }
          MPI_Bcast(&key_size, 1, MPI_UNSIGNED, 0, communicator);
          if (this_mpi_process != 0)
            {
              keys[i].resize(key_size);
            }
          MPI_Bcast(const_cast<char *>(keys[i].data()),
                    key_size,
                    MPI_CHAR,
                    0,
                    communicator);
        }

      // Number of particles in the file
      unsigned int n_particles_in_file = 0;
      if (this_mpi_process == 0)
        {
          n_particles_in_file = particles_data.at(keys[0]).size();
        }
      MPI_Bcast(&n_particles_in_file, 1, MPI_UNSIGNED, 0, communicator);

      // Limit the number of particles to insert
      unsigned int n_particles_to_insert =
        std::min(remaining_particles_of_each_type, n_particles_in_file);

      // Distribute particles among ranks
      const int          n_mpi_procs = Utilities::MPI::n_mpi_processes(communicator);
      const unsigned int n_particles_per_proc = n_particles_to_insert / n_mpi_procs;
      const unsigned int n_excess_particles   = n_particles_to_insert % n_mpi_procs;

      unsigned int my_n_particles = n_particles_per_proc;
      unsigned int my_first_p     = this_mpi_process * n_particles_per_proc;

      if ((unsigned int)this_mpi_process < n_excess_particles)
        {
          my_n_particles++;
          my_first_p += this_mpi_process;
        }
      else
        {
          my_first_p += n_excess_particles;
        }

      // Prepare local particles data
      std::map<std::string, std::vector<double>> local_particles_data;
      for (const auto &key : keys)
        {
          local_particles_data[key].resize(my_n_particles);
          // Scatter the data for this key
          std::vector<int>  sendcounts(n_mpi_procs);
          std::vector<int>  displs(n_mpi_procs);
          const double     *sendbuf = nullptr;
          if (this_mpi_process == 0)
            {
              sendbuf = particles_data.at(key).data();
              for (int i = 0; i < n_mpi_procs; ++i)
                {
                  sendcounts[i] = n_particles_per_proc;
                  displs[i]     = i * n_particles_per_proc;
                  if ((unsigned int)i < n_excess_particles)
                    {
                      sendcounts[i]++;
                      displs[i] += i;
                    }
                  else
                    {
                      displs[i] += n_excess_particles;
                    }
                }
            }
          MPI_Scatterv(sendbuf,
                       sendcounts.data(),
                       displs.data(),
                       MPI_DOUBLE,
                       local_particles_data[key].data(),
                       my_n_particles,
                       MPI_DOUBLE,
                       0,
                       communicator);
        }

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      std::vector<Point<dim>>          all_insertion_points_on_proc;
      std::vector<std::vector<double>> all_particle_properties;

      for (unsigned int p = 0; p < my_n_particles; ++p)
        {
          if constexpr (dim == 2)
            {
              all_insertion_points_on_proc.emplace_back(Point<dim>(
                {local_particles_data["p_x"][p], local_particles_data["p_y"][p]}));
            }

          if constexpr (dim == 3)
            {
              all_insertion_points_on_proc.emplace_back(Point<dim>(
                {local_particles_data["p_x"][p],
                 local_particles_data["p_y"][p],
                 local_particles_data["p_z"][p]}));
            }
        }

      // Assign inserted particles properties
      this->assign_particle_properties_for_file_insertion(
        dem_parameters,
        my_n_particles,
        0,
        local_particles_data,
        all_particle_properties);

      // Insert all particles at once to avoid communicator mismatches
      particle_handler.insert_global_particles(all_insertion_points_on_proc,
                                               global_bounding_boxes,
                                               all_particle_properties);

      remaining_particles_of_each_type -= n_particles_to_insert;

      ConditionalOStream pcout(std::cout, this_mpi_process == 0);
      this->print_insertion_info(n_particles_to_insert,
                                 remaining_particles_of_each_type,
                                 this->current_inserting_particle_type,
                                 pcout);
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionFile<dim, PropertiesIndex>::
  assign_particle_properties_for_file_insertion(
    const DEMSolverParameters<dim>             &dem_parameters,
    const unsigned int                         &inserted_this_step_this_proc,
    const unsigned int                         &start_index,
    std::map<std::string, std::vector<double>> &particles_data,
    std::vector<std::vector<double>>           &particle_properties)
{
  // Clearing and resizing particle_properties
  particle_properties.reserve(inserted_this_step_this_proc);

  // Getting properties as local parameters
  auto physical_properties = dem_parameters.lagrangian_physical_properties;

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0;
       particle_counter < inserted_this_step_this_proc;
       ++particle_counter)
    {
      unsigned int global_particle_counter = start_index + particle_counter;
      double       type = this->current_inserting_particle_type;
      double       diameter =
        particles_data["diameters"][global_particle_counter];
      double density =
        physical_properties
          .density_particle[this->current_inserting_particle_type];
      double vel_x   = particles_data["v_x"][global_particle_counter];
      double vel_y   = particles_data["v_y"][global_particle_counter];
      double vel_z   = particles_data["v_z"][global_particle_counter];
      double omega_x = particles_data["w_x"][global_particle_counter];
      double omega_y = particles_data["w_y"][global_particle_counter];
      double omega_z = particles_data["w_z"][global_particle_counter];
      double mass    = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);

      std::vector<double> properties_of_one_particle{
        type, diameter, mass, vel_x, vel_y, vel_z, omega_x, omega_y, omega_z};

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::DEMMPProperties::PropertiesIndex>)
        {
          double T = particles_data["T"][global_particle_counter];
          double specific_heat =
            physical_properties
              .specific_heat_particle[this->current_inserting_particle_type];
          properties_of_one_particle.push_back(T);
          properties_of_one_particle.push_back(specific_heat);
        }

      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::CFDDEMProperties::PropertiesIndex>)
        {
          // Push back all zero variables for the CFD-DEM coupling properties
          for (unsigned int i = properties_of_one_particle.size();
               i < PropertiesIndex::n_properties;
               ++i)
            properties_of_one_particle.push_back(0.);
        }

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}

template class InsertionFile<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionFile<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionFile<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::DEMMPProperties::PropertiesIndex>;
