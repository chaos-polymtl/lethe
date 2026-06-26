// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/utilities.h>

#include <dem/insertion_volume.h>

using namespace DEM;

// The constructor of the volume insertion class. In the constructor, we
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
  , acceptance_fct(
      dem_parameters.insertion_information.insertion_acceptance_fct)
{
  // Initializing the current inserting particle type
  current_inserting_particle_type = 0;
  this->maximum_diameter          = maximum_particle_diameter;

  // We need to fill the map filted_id_to_box_id.
  // To do this, we loop of every insertion point inside the box considering
  // that the acceptance_fct accepts every point.
  // For each points, we check if it respects the condition. If so, we insert
  // the ID associated with this point inside the map.
  set_filtered_id_map(dem_parameters.insertion_info);
}

// The main insertion function. Insert_global_function is used to insert
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

  // Check to see if the remaining uninserted particles are equal to zero or
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

      unsigned int inserted_this_step;
      this->calculate_insertion_domain_maximum_particle_number(
        dem_parameters.insertion_info, pcout, inserted_this_step);

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // Distributing particles between processors
      unsigned int inserted_this_step_this_proc =
        static_cast<unsigned int>(floor(inserted_this_step / n_mpi_process));
      if (this_mpi_process == (n_mpi_process - 1))
        inserted_this_step_this_proc = static_cast<unsigned int>(
          inserted_this_step -
          (n_mpi_process - 1) * floor(inserted_this_step / n_mpi_process));

      // Call the random number generator
      std::vector<double> random_number_vector;
      random_number_vector.reserve(inserted_this_step_this_proc);
      create_random_number_container(
        random_number_vector,
        inserted_this_step,
        dem_parameters.insertion_information.insertion_maximum_offset,
        dem_parameters.insertion_information.seed_for_insertion);

      Point<dim>              insertion_location;
      std::vector<Point<dim>> insertion_points_on_proc;
      insertion_points_on_proc.reserve(inserted_this_step_this_proc);

      // Find the first and the last particle id for each process
      // The number of particles on the last process is different
      unsigned int first_id;
      unsigned int last_id;
      if (this_mpi_process == (n_mpi_process - 1))
        {
          first_id = inserted_this_step - inserted_this_step_this_proc;
          last_id  = inserted_this_step;
        }
      // For the processes 1 : n-1
      else
        {
          first_id = this_mpi_process * inserted_this_step_this_proc;
          last_id  = (this_mpi_process + 1) * inserted_this_step_this_proc;
        }

      // Looping through the particles on each process and finding their
      // insertion location
      unsigned int particle_counter = 0;
      for (unsigned int id = first_id; id < last_id; ++id, ++particle_counter)
        {
          find_insertion_location(
            insertion_location,
            filted_id_to_box_id.at(id),
            random_number_vector[particle_counter],
            random_number_vector[inserted_this_step - particle_counter - 1],
            dem_parameters.insertion_information);

          insertion_points_on_proc.push_back(insertion_location);
        }

      std::vector<std::vector<double>> particle_properties;

      // Assigning inserted particles properties using
      // assign_particle_properties function
      this->assign_particle_properties(dem_parameters,
                                       inserted_this_step_this_proc,
                                       current_inserting_particle_type,
                                       insertion_points_on_proc,
                                       particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
                                               global_bounding_boxes,
                                               particle_properties);

      // Updating remaining particles
      particles_of_each_type_remaining -= inserted_this_step;

      this->print_insertion_information(inserted_this_step,
                                        particles_of_each_type_remaining,
                                        current_inserting_particle_type,
                                        pcout);
    }
}

// This function assigns the insertion points of the inserted particles
template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::find_insertion_location(
  Point<dim>        &insertion_location,
  const unsigned int id,
  const double       random_number1,
  const double       random_number2,
  const Parameters::Lagrangian::InsertionInfo<dim>
    &insertion_informationrmation)
{
  std::vector<int> insertion_index;
  insertion_index.resize(dim);

  unsigned int axis_0, axis_1;
  int          number_of_particles_0, number_of_particles_1;

  // First direction (axis) to have particles inserted
  axis_0                = insertion_informationrmation.direction_sequence.at(0);
  number_of_particles_0 = this->number_of_particles_directions[axis_0];
  insertion_index[axis_0] = id % number_of_particles_0;
  insertion_location[axis_0] =
    this->axis_min[axis_0] +
    ((insertion_index[axis_0] + 0.5) *
       insertion_informationrmation.distance_threshold -
     random_number1) *
      this->maximum_diameter;

  // Second direction (axis) to have particles inserted
  axis_1 = insertion_informationrmation.direction_sequence.at(1);

  if constexpr (dim == 2)
    {
      insertion_index[axis_1] = static_cast<int>(id / number_of_particles_0);
      insertion_location[axis_1] =
        this->axis_min[axis_1] +
        ((insertion_index[axis_1] + 0.5) *
           insertion_informationrmation.distance_threshold -
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
        this->axis_min[axis_1] +
        ((insertion_index[axis_1] + 0.5) *
           insertion_informationrmation.distance_threshold -
         random_number2) *
          this->maximum_diameter;

      // Third direction (axis) to have particles inserted
      unsigned int axis_2;
      axis_2 = insertion_informationrmation.direction_sequence.at(2);
      insertion_index[axis_2] =
        static_cast<int>(id / (number_of_particles_0 * number_of_particles_1));
      insertion_location[axis_2] =
        this->axis_min[axis_2] +
        ((insertion_index[axis_2] + 0.5) *
           insertion_informationrmation.distance_threshold -
         random_number1) *
          this->maximum_diameter;
    }
}


template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::set_filtered_id_map(
  const InsertionInfo<dim> &insertion_informationrmation)
{
  // Checking if the insertion directions are valid (no repetition)
  int axis_sum = 0;
  if constexpr (dim == 2)
    {
      axis_sum = insertion_informationrmation.direction_sequence[0] +
                 insertion_informationrmation.direction_sequence[1];

      AssertThrow(
        axis_sum == 1,
        ExcMessage("Invalid insertion directions: 2 directions are the same "));
    }
  if constexpr (dim == 3)
    {
      axis_sum = insertion_informationrmation.direction_sequence[0] +
                 insertion_informationrmation.direction_sequence[1] +
                 insertion_informationrmation.direction_sequence[2];

      AssertThrow(
        axis_sum == 3,
        ExcMessage(
          "Invalid insertion directions: at least 2 directions are the same "));
    }

  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box before the acceptance function.
  int maximum_particle_number = 1;

  number_of_particles_directions.resize(dim);
  axis_min.resize(dim);
  axis_max.resize(dim);

  std::vector<unsigned int> axis_list = {
    insertion_informationrmation.direction_sequence[0],
    insertion_informationrmation.direction_sequence[1]};

  if constexpr (dim == 3)
    {
      axis_list.push_back(insertion_informationrmation.direction_sequence[2]);
    }

  // Assigning the minimum and maximum positions of the insertion box in respect
  // to the axis order
  for (unsigned int axis : axis_list)
    {
      switch (axis)
        {
          case 0:
            axis_min[0] = insertion_informationrmation.insertion_box_point_1(0);
            axis_max[0] = insertion_informationrmation.insertion_box_point_2(0);
            break;
          case 1:
            axis_min[1] = insertion_informationrmation.insertion_box_point_1(1);
            axis_max[1] = insertion_informationrmation.insertion_box_point_2(1);
            break;
          case 2:
            axis_min[2] = insertion_informationrmation.insertion_box_point_1(2);
            axis_max[2] = insertion_informationrmation.insertion_box_point_2(2);
            break;
          default:
            AssertThrow(false,
                        ExcMessage("Insertion direction must be 0, 1 or 2"));
        }

      // Assign max number of particles according to the direction and calculate
      // the total max number (maximum_particle_number = max_x * max_y * max_z)
      int number_of_particles =
        static_cast<int>((axis_max[axis] - axis_min[axis]) /
                         (insertion_informationrmation.distance_threshold *
                          this->maximum_diameter));
      number_of_particles_directions[axis] = number_of_particles;

      maximum_particle_number *= number_of_particles;
    }
  maximum_particle_number =
    std::min(maximum_particle_number,
             dem_parameters.insertion_information.inserted_this_step);
  // Now, we know that the ID before the acceptance function will go from 0
  // to maximum_particle_number - 1 .
  // We count the number of insertion points that respect the acceptance
  // function and store the ID of those points

  unsigned int filtered_id_count = 0;
  Point<dim>   insertion_location;
  for (int id = 0; id < maximum_particle_number; ++id)
    {
      // Create the insertion point
      find_insertion_location(
        insertion_location, id, 0., 0., insertion_informationrmation);

      if (acceptance_fct->value(insertion_location) > 0.)
        {
          filted_id_to_box_id.insert(std::make_pair(filtered_id_count, id));
          ++filtered_id_count;
        }
    }
}

template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::
  calculate_insertion_domain_maximum_particle_number(
    const InsertionInfo<dim> &insertion_information,
    const ConditionalOStream &pcout,
    unsigned int             &inserted_this_step)
{
  // Maximum number of particles that fit inside the insertion box filter by
  // the function.
  unsigned int maximum_particle_number = filted_id_to_box_id.size();

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number)
    {
      pcout << "Warning: the requested number of particles for insertion ("
            << insertion_information.inserted_this_step
            << ") is higher than maximum expected number of particles ("
            << maximum_particle_number << ")" << std::endl;

      // Updating the number of inserted particles at each step
      inserted_this_step = maximum_particle_number;
    }
  else
    {
      inserted_this_step = insertion_information.inserted_this_step;
    }

  // The inserted_this_step value is the minimum of
  // particles_of_each_type_remaining and inserted_this_step
  inserted_this_step =
    std::min(particles_of_each_type_remaining, inserted_this_step);
}

template class InsertionVolume<2, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<2, DEM::DEMMPProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::CFDDEMProperties::PropertiesIndex>;
template class InsertionVolume<3, DEM::DEMMPProperties::PropertiesIndex>;
