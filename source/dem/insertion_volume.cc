// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

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
      this->inserted_this_step_this_proc = static_cast<unsigned int>(
        floor(this->inserted_this_step / n_mpi_process));
      if (this_mpi_process == (n_mpi_process - 1))
        this->inserted_this_step_this_proc = static_cast<unsigned int>(
          this->inserted_this_step -
          (n_mpi_process - 1) *
            floor(this->inserted_this_step / n_mpi_process));

      // Call random number generator
      std::vector<double> random_number_vector;
      random_number_vector.reserve(this->inserted_this_step_this_proc);
      create_random_number_container(
        random_number_vector,
        this->inserted_this_step,
        dem_parameters.insertion_info.insertion_maximum_offset,
        dem_parameters.insertion_info.seed_for_insertion);

      Point<dim>              insertion_locations;
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
          find_insertion_location(
            insertion_locations,
            id,
            random_number_vector[particle_counter],
            random_number_vector[this->inserted_this_step - particle_counter -
                                 1],
            dem_parameters.insertion_info);

          insertion_points_on_proc.push_back(insertion_locations);
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

template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::find_insertion_location(
  Point<dim>                                       &insertion_location,
  const unsigned int                                id,
  const double                                      random_number1,
  const double                                      random_number2,
  const Parameters::Lagrangian::InsertionInfo<dim> &insertion_information)
{
  // Map random numbers to their corresponding insertion sequence order
  // axis_sequence[0] uses random_number1, axis_sequence[1] uses random_number2,
  // axis_sequence[2] uses random_number1
  const double random_numbers[3] = {random_number1,
                                    random_number2,
                                    random_number1};

  // Unflatten the 1D 'id' into structured grid indices based on the direction
  // sequence
  std::vector<int> insertion_index(dim);
  unsigned int     remainder = id;

  for (int i = 0; i < dim; ++i)
    {
      const unsigned int axis = insertion_information.direction_sequence.at(i);
      const int particles_in_axis = this->number_of_particles_directions[axis];

      // Extract the coordinate index for this specific axis
      if (i < dim - 1)
        {
          insertion_index[axis] = remainder % particles_in_axis;
          remainder /= particles_in_axis;
        }
      else
        {
          // For the final axis, the remaining quotient is the index
          insertion_index[axis] = remainder;
        }

      // Calculate the spatial location for this axis
      insertion_location[axis] =
        this->axis_min[axis] + ((insertion_index[axis] + 0.5) *
                                  insertion_information.distance_threshold -
                                random_numbers[i]) *
                                 this->maximum_diameter;
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::
  calculate_insertion_domain_maximum_particle_number(
    const DEMSolverParameters<dim> &dem_parameters,
    const ConditionalOStream       &pcout)
{
  // Getting properties as local parameters
  const auto insertion_information = dem_parameters.insertion_info;

  // Checking if the insertion directions are valid (no repetition)
  int axis_sum = 0;
  if constexpr (dim == 2)
    {
      axis_sum = insertion_information.direction_sequence[0] +
                 insertion_information.direction_sequence[1];

      AssertThrow(
        axis_sum == 1,
        ExcMessage("Invalid insertion directions: 2 directions are the same "));
    }
  if constexpr (dim == 3)
    {
      axis_sum = insertion_information.direction_sequence[0] +
                 insertion_information.direction_sequence[1] +
                 insertion_information.direction_sequence[2];

      AssertThrow(
        axis_sum == 3,
        ExcMessage(
          "Invalid insertion directions: at least 2 directions are the same "));
    }


  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number = 1;

  number_of_particles_directions.resize(dim);
  axis_min.resize(dim);
  axis_max.resize(dim);

  std::vector<unsigned int> axis_list = {
    insertion_information.direction_sequence[0],
    insertion_information.direction_sequence[1]};

  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_information.direction_sequence[2]);
    }

  // Assigning the minimum and maximum positions of the insertion box in respect
  // to the axis order
  for (unsigned int axis : axis_list)
    {
      switch (axis)
        {
          case 0:
            axis_min[0] = insertion_information.insertion_box_point_1(0);
            axis_max[0] = insertion_information.insertion_box_point_2(0);
            break;
          case 1:
            axis_min[1] = insertion_information.insertion_box_point_1(1);
            axis_max[1] = insertion_information.insertion_box_point_2(1);
            break;
          case 2:
            axis_min[2] = insertion_information.insertion_box_point_1(2);
            axis_max[2] = insertion_information.insertion_box_point_2(2);
            break;
          default:
            AssertThrow(false,
                        ExcMessage("Insertion direction must be 0, 1 or 2"));
        }

      // Assign max number of particles according to the direction and calculate
      // the total max number (maximum_particle_number = max_x * max_y * max_z)
      int number_of_particles = static_cast<int>(
        (axis_max[axis] - axis_min[axis]) /
        (insertion_information.distance_threshold * this->maximum_diameter));
      number_of_particles_directions[axis] = number_of_particles;
      maximum_particle_number *= number_of_particles;
    }

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number)
    {
      pcout << "Warning: the requested number of particles for insertion ("
            << insertion_information.inserted_this_step
            << ") is higher than maximum expected number of particles ("
            << maximum_particle_number << ")" << std::endl;

      // Updating the number of inserted particles at each step
      this->inserted_this_step = maximum_particle_number;
    }
  else
    {
      this->inserted_this_step = insertion_information.inserted_this_step;
    }
}

template class InsertionVolume<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMMPProperties::PropertiesIndex>;
