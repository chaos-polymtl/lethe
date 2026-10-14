// SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/dem_properties.h>
#include <core/utilities.h>
#include <core/vtu_reader.h>

#include <dem/insertion_file.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/mpi.h>

#include <filesystem>

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

      MPI_Comm communicator = triangulation.get_mpi_communicator();
      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

      // Read the input file. The particles of a .vtu or .pvtu file are read by
      // all processes, each reading part of the .vtu files, while the
      // particles of a text file are only inserted by processor 0.
      std::map<std::string, std::vector<double>> particles_data;
      const std::string &file_name = insertion_files.at(current_file_id);
      const std::string  file_extension =
        std::filesystem::path(file_name).extension().string();
      unsigned int n_particles_in_file_this_proc;
      if (file_extension == ".vtu" || file_extension == ".pvtu")
        {
          read_vtu_particles_data(file_name, communicator, particles_data);
          n_particles_in_file_this_proc = particles_data["p_x"].size();
        }
      else
        {
          fill_vectors_from_file(particles_data, file_name, ";");
          n_particles_in_file_this_proc =
            this_mpi_process == 0 ? particles_data["p_x"].size() : 0;
        }
      current_file_id++;
      current_file_id = current_file_id % number_of_files;

      // Number of particles in the file read by the processors of lower rank
      // and number of particles in the file
      const auto [n_particles_in_file_lower_procs, n_particles_in_file] =
        Utilities::MPI::partial_and_total_sum(n_particles_in_file_this_proc,
                                              communicator);

      // Adjusting the value in case we exceed the maximum number of particle in
      // the simulation. Only the first particles of the file are inserted.
      const unsigned int n_total_particles_to_insert =
        std::min(remaining_particles_of_each_type, n_particles_in_file);
      const unsigned int n_particles_to_insert_this_proc =
        std::min(n_particles_in_file_this_proc,
                 n_total_particles_to_insert -
                   std::min(n_particles_in_file_lower_procs,
                            n_total_particles_to_insert));

      std::vector<Point<dim>> insertion_points_on_proc_this_step;
      insertion_points_on_proc_this_step.reserve(
        n_particles_to_insert_this_proc);

      for (unsigned int p = 0; p < n_particles_to_insert_this_proc; ++p)
        {
          if constexpr (dim == 2)
            {
              insertion_points_on_proc_this_step.emplace_back(Point<dim>(
                {particles_data["p_x"][p], particles_data["p_y"][p]}));
            }

          if constexpr (dim == 3)
            {
              insertion_points_on_proc_this_step.emplace_back(
                Point<dim>({particles_data["p_x"][p],
                            particles_data["p_y"][p],
                            particles_data["p_z"][p]}));
            }
        }

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // A vector of vectors, which contains all the properties of all particles
      // about to get inserted
      std::vector<std::vector<double>> particle_properties;

      // Assign inserted particles properties
      this->assign_particle_properties_for_file_insertion(
        dem_parameters,
        n_particles_to_insert_this_proc,
        particles_data,
        particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(
        insertion_points_on_proc_this_step,
        global_bounding_boxes,
        particle_properties);

      // Update number of particle remaining to be inserted
      remaining_particles_of_each_type -= n_total_particles_to_insert;


      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
      this->print_insertion_info(n_total_particles_to_insert,
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
      double type     = this->current_inserting_particle_type;
      double diameter = particles_data["diameters"][particle_counter];
      double density =
        physical_properties
          .density_particle[this->current_inserting_particle_type];
      double vel_x   = particles_data["v_x"][particle_counter];
      double vel_y   = particles_data["v_y"][particle_counter];
      double vel_z   = particles_data["v_z"][particle_counter];
      double omega_x = particles_data["w_x"][particle_counter];
      double omega_y = particles_data["w_y"][particle_counter];
      double omega_z = particles_data["w_z"][particle_counter];
      double mass    = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);

      // Every property which is not set explicitly below, such as the CFD-DEM
      // coupling properties, is left to zero.
      std::vector<double> properties_of_one_particle(
        PropertiesIndex::n_properties, 0.);

      properties_of_one_particle[PropertiesIndex::type]    = type;
      properties_of_one_particle[PropertiesIndex::dp]      = diameter;
      properties_of_one_particle[PropertiesIndex::mass]    = mass;
      properties_of_one_particle[PropertiesIndex::v_x]     = vel_x;
      properties_of_one_particle[PropertiesIndex::v_y]     = vel_y;
      properties_of_one_particle[PropertiesIndex::v_z]     = vel_z;
      properties_of_one_particle[PropertiesIndex::omega_x] = omega_x;
      properties_of_one_particle[PropertiesIndex::omega_y] = omega_y;
      properties_of_one_particle[PropertiesIndex::omega_z] = omega_z;

      if constexpr (DEM::has_thermal_properties<PropertiesIndex>)
        {
          properties_of_one_particle[PropertiesIndex::T] =
            particles_data["T"][particle_counter];
          properties_of_one_particle[PropertiesIndex::specific_heat] =
            physical_properties
              .specific_heat_particle[this->current_inserting_particle_type];
        }

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionFile<dim, PropertiesIndex>::read_vtu_particles_data(
  const std::string                          &file_name,
  const MPI_Comm                             &mpi_communicator,
  std::map<std::string, std::vector<double>> &particles_data)
{
  // Arrays of the .vtu files to read and keys under which their components
  // are stored. The names of the arrays are the ones used to write the
  // particles, and the keys are the column names of the insertion text files.
  const std::vector<std::pair<std::string, int>> properties_name =
    ParticleProperties<dim, PropertiesIndex>::get_properties_name();
  const std::vector<std::string> coordinate_keys = {"p_x", "p_y", "p_z"};

  std::map<std::string, std::vector<std::string>> array_keys;
  array_keys["Points"] =
    std::vector<std::string>(coordinate_keys.begin(),
                             coordinate_keys.begin() + dim);
  array_keys[properties_name[PropertiesIndex::v_x].first]     = {"v_x",
                                                                 "v_y",
                                                                 "v_z"};
  array_keys[properties_name[PropertiesIndex::omega_x].first] = {"w_x",
                                                                 "w_y",
                                                                 "w_z"};
  array_keys[properties_name[PropertiesIndex::dp].first]      = {"diameters"};
  if constexpr (DEM::has_thermal_properties<PropertiesIndex>)
    array_keys[properties_name[PropertiesIndex::T].first] = {"T"};

  // Every process lists the .vtu files, since the .pvtu file is small
  const std::vector<std::string> vtu_file_names = get_vtu_file_names(file_name);

  // Each process reads a contiguous block of .vtu files
  const IndexSet vtu_files_this_proc =
    Utilities::MPI::create_evenly_distributed_partitioning(
      mpi_communicator, vtu_file_names.size());
  for (const auto vtu_file_id : vtu_files_this_proc)
    fill_vectors_from_vtu(particles_data,
                          vtu_file_names[vtu_file_id],
                          array_keys);
}

template class InsertionFile<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionFile<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionFile<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionFile<2, DEM::CFDDEMMPProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionFile<3, DEM::CFDDEMMPProperties::PropertiesIndex>;
