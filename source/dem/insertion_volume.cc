// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <dem/insertion_volume.h>

using namespace DEM;

// The constructor of volume insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim, typename PropertiesIndex>
InsertionVolume<dim, PropertiesIndex>::InsertionVolume(
  const std::vector<std::shared_ptr<Distribution>>
    &size_distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters,
  const double                                     maximum_particle_diameter)
  : Insertion<dim, PropertiesIndex>(size_distribution_object_container,
                                    triangulation,
                                    dem_parameters)
  , particles_of_each_type_remaining(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;
  this->inserted_this_step        = 0;
  this->maximum_diameter          = maximum_particle_diameter;
}

// The main insertion function. Insert_global_function is utilized to insert
// the particles
template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  if (particles_of_each_type_remaining == 0 &&
      current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      particles_of_each_type_remaining =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }

  // Check to see if the remaining un-inserted particles are equal to zero or
  // not
  if (particles_of_each_type_remaining != 0)
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

      MPI_Comm           communicator = triangulation.get_mpi_communicator();
      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      auto n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);


      this->calculate_insertion_domain_maximum_particle_number(dem_parameters,
                                                               pcout);
      // The inserted_this_step value is the minimum of
      // particles_of_each_type_remaining and inserted_this_step
      this->inserted_this_step =
        std::min(particles_of_each_type_remaining, this->inserted_this_step);

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // Distributing particles between processors
      this->inserted_this_step_this_proc =
        floor(this->inserted_this_step / n_mpi_process);
      if (this_mpi_process == (n_mpi_process - 1))
        this->inserted_this_step_this_proc =
          this->inserted_this_step -
          (n_mpi_process - 1) * floor(this->inserted_this_step / n_mpi_process);

      // Call random number generator
      std::vector<double> random_number_vector;
      random_number_vector.reserve(this->inserted_this_step_this_proc);
      this->create_random_number_container(
        random_number_vector,
        dem_parameters.insertion_info.insertion_maximum_offset,
        dem_parameters.insertion_info.seed_for_insertion);

      Point<dim>              insertion_location;
      std::vector<Point<dim>> insertion_points_on_proc;
      insertion_points_on_proc.reserve(this->inserted_this_step_this_proc);

      // Find the first and the last particle id for each process
      // The number of particles on the last process is different
      unsigned int first_id;
      unsigned int last_id;
      if (this_mpi_process == (n_mpi_process - 1))
        {
          first_id =
            this->inserted_this_step - this->inserted_this_step_this_proc;
          last_id = this->inserted_this_step;
        }
      // For the processes 1 : n-1
      else
        {
          first_id = this_mpi_process * this->inserted_this_step_this_proc;
          last_id = (this_mpi_process + 1) * this->inserted_this_step_this_proc;
        }

      // Looping through the particles on each process and finding their
      // insertion location
      unsigned int particle_counter = 0;
      for (unsigned int id = first_id; id < last_id; ++id, ++particle_counter)
        {
          find_insertion_location_volume(
            insertion_location,
            id,
            random_number_vector[particle_counter],
            random_number_vector[this->inserted_this_step - particle_counter -
                                 1],
            dem_parameters.insertion_info);
          insertion_points_on_proc.push_back(insertion_location);
        }

      std::vector<std::vector<double>> particle_properties;

      // Assigning inserted particles properties using
      // assign_particle_properties function
      this->assign_particle_properties(dem_parameters,
                                       this->inserted_this_step_this_proc,
                                       current_inserting_particle_type,
                                       insertion_points_on_proc,
                                       particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
                                               global_bounding_boxes,
                                               particle_properties);

      // Updating remaining particles
      particles_of_each_type_remaining -= this->inserted_this_step;

      this->print_insertion_info(this->inserted_this_step,
                                 particles_of_each_type_remaining,
                                 current_inserting_particle_type,
                                 pcout);
    }
}

// This function creates a vector of random doubles using the input parameters
// in the parameter handler
template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::create_random_number_container(
  std::vector<double> &random_container,
  const double         maximum_range,
  const int            seed_for_insertion)
{
  for (unsigned int i = 0; i < this->inserted_this_step; ++i)
    {
      srand(seed_for_insertion * (i + 1));
      random_container.push_back((((double)rand()) / ((double)RAND_MAX)) *
                                 maximum_range);
    }
}

// This function assigns the insertion points of the inserted particles
template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::find_insertion_location_volume(
  Point<dim>                                       &insertion_location,
  const unsigned int                                id,
  const double                                      random_number1,
  const double                                      random_number2,
  const Parameters::Lagrangian::InsertionInfo<dim> &insertion_information)
{
  std::vector<int> insertion_index;
  insertion_index.resize(dim);

  unsigned int axis_0, axis_1;
  int          number_of_particles_0, number_of_particles_1;

  // First direction (axis) to have particles inserted
  axis_0                  = insertion_information.direction_sequence.at(0);
  number_of_particles_0   = this->number_of_particles_directions[axis_0];
  insertion_index[axis_0] = id % number_of_particles_0;
  insertion_location[axis_0] =
    this->axis_min[axis_0] + ((insertion_index[axis_0] + 0.5) *
                                insertion_information.distance_threshold -
                              random_number1) *
                               this->maximum_diameter;

  // Second direction (axis) to have particles inserted
  axis_1 = insertion_information.direction_sequence.at(1);

  if constexpr (dim == 2)
    {
      insertion_index[axis_1] = static_cast<int>(id / number_of_particles_0);
      insertion_location[axis_1] =
        this->axis_min[axis_1] + ((insertion_index[axis_1] + 0.5) *
                                    insertion_information.distance_threshold -
                                  random_number2) *
                                   this->maximum_diameter;
    }
  else
    {
      number_of_particles_1 = this->number_of_particles_directions[axis_1];
      insertion_index[axis_1] =
        static_cast<int>(id % (number_of_particles_0 * number_of_particles_1)) /
        (number_of_particles_0);
      insertion_location[axis_1] =
        this->axis_min[axis_1] + ((insertion_index[axis_1] + 0.5) *
                                    insertion_information.distance_threshold -
                                  random_number2) *
                                   this->maximum_diameter;

      // Third direction (axis) to have particles inserted
      unsigned int axis_2;
      axis_2 = insertion_information.direction_sequence.at(2);
      insertion_index[axis_2] =
        static_cast<int>(id / (number_of_particles_0 * number_of_particles_1));
      insertion_location[axis_2] =
        this->axis_min[axis_2] + ((insertion_index[axis_2] + 0.5) *
                                    insertion_information.distance_threshold -
                                  random_number1) *
                                   this->maximum_diameter;
    }
}

template class InsertionVolume<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMMPProperties::PropertiesIndex>;
