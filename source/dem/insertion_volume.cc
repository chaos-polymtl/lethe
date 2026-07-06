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
  , acceptance_fct(dem_parameters.insertion_info.insertion_acceptance_fct)
  , maximum_diameter(maximum_particle_diameter)
{
  // Initializing the current inserting particle type
  current_inserting_particle_type = 0;

  // We need to fill the vector filted_box_index.
  // To do this, we loop of every insertion point inside the box considering
  // that the acceptance_fct accepts every point.
  // For each point, we check if it respects the condition. If so, we insert
  // the index associated with this point inside the vector.
  set_filtered_index(dem_parameters.insertion_info);
}


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

  // Check if the remaining number of uninserted particles is equal to zero
  if (particles_of_each_type_remaining != 0)
    {
      // Remove the particle if the remove box region feature is activated
      if (this->removing_particles_in_region)
        {
          if (this->mark_for_update)
            {
              this->find_cells_in_removing_box(triangulation);
              this->mark_for_update = false;
            }
          this->remove_particles_in_box(particle_handler);
        }

      // Get the MPI communicator and the parallel cout
      MPI_Comm           communicator = triangulation.get_mpi_communicator();
      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      // Compute the number of particle to insert during this step.
      unsigned int inserted_this_step;
      adjust_insertion_count_by_insertion_box_capacity(
        dem_parameters.insertion_info, pcout, inserted_this_step);

      // Get global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // Call the random number generator for the offsets
      std::vector<double> random_number_vector;
      random_number_vector.reserve(filted_box_index.size());
      create_random_number_container(
        random_number_vector,
        filted_box_index.size(),
        dem_parameters.insertion_info.insertion_maximum_offset,
        dem_parameters.insertion_info.seed_for_insertion);

      Point<dim>              insertion_location;
      std::vector<Point<dim>> insertion_points_on_proc;
      insertion_points_on_proc.reserve(filted_box_index.size());

      // Looping through the particles on each process and finding their
      // insertion location
      unsigned int particle_counter = 0;
      for (unsigned int global_index = first_index_this_proc;
           global_index < last_index_this_proc &&
           global_index < inserted_this_step;
           ++global_index, ++particle_counter)
        {
          find_insertion_location(insertion_location,
                                  filted_box_index.at(particle_counter),
                                  random_number_vector.at(particle_counter),
                                  random_number_vector.at(
                                    filted_box_index.size() - particle_counter -
                                    1),
                                  dem_parameters.insertion_info);

          insertion_points_on_proc.push_back(insertion_location);
        }

      std::vector<std::vector<double>> particle_properties;

      // Assigning inserted particles properties using
      // assign_particle_properties function
      this->assign_particle_properties(dem_parameters,
                                       particle_counter,
                                       current_inserting_particle_type,
                                       insertion_points_on_proc,
                                       particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
                                               global_bounding_boxes,
                                               particle_properties);

      // Updating remaining particles
      particles_of_each_type_remaining -= inserted_this_step;

      this->print_insertion_info(inserted_this_step,
                                 particles_of_each_type_remaining,
                                 current_inserting_particle_type,
                                 pcout);
    }
}

// This function assigns the insertion points of the inserted particles
template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::find_insertion_location(
  Point<dim>                                       &insertion_location,
  const unsigned int                                id,
  const double                                      random_number1,
  const double                                      random_number2,
  const Parameters::Lagrangian::InsertionInfo<dim> &insertion_information)
{
  std::vector<int> insertion_index;
  insertion_index.resize(dim);

  // First direction (axis) to have particles inserted
  unsigned int axis_0 = insertion_information.direction_sequence[0];
  unsigned int number_of_particles_0 = number_of_particles_directions[axis_0];
  insertion_index[axis_0]            = id % number_of_particles_0;
  insertion_location[axis_0] =
    axis_min[axis_0] + ((insertion_index[axis_0] + 0.5) *
                          insertion_information.distance_threshold -
                        random_number1) *
                         maximum_diameter;

  // Second direction (axis) to have particles inserted
  unsigned int axis_1 = insertion_information.direction_sequence.at(1);

  if constexpr (dim == 2)
    {
      insertion_index[axis_1] = static_cast<int>(id / number_of_particles_0);
      insertion_location[axis_1] =
        axis_min[axis_1] + ((insertion_index[axis_1] + 0.5) *
                              insertion_information.distance_threshold -
                            random_number2) *
                             maximum_diameter;
    }
  else
    {
      unsigned int number_of_particles_1 =
        number_of_particles_directions[axis_1];
      insertion_index[axis_1] =
        static_cast<int>(id % (number_of_particles_0 * number_of_particles_1)) /
        (number_of_particles_0);
      insertion_location[axis_1] =
        axis_min[axis_1] + ((insertion_index[axis_1] + 0.5) *
                              insertion_information.distance_threshold -
                            random_number2) *
                             maximum_diameter;

      // Third direction (axis) to have particles inserted
      unsigned int axis_2 = insertion_information.direction_sequence[2];
      insertion_index[axis_2] =
        static_cast<int>(id / (number_of_particles_0 * number_of_particles_1));
      insertion_location[axis_2] =
        axis_min[axis_2] + ((insertion_index[axis_2] + 0.5) *
                              insertion_information.distance_threshold -
                            random_number1) *
                             maximum_diameter;
    }
}


template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::set_filtered_index(
  const InsertionInfo<dim> &insertion_information)
{
  // If you don't insert particles in this simulation, we don't need to build
  // the container.
  if (insertion_information.inserted_this_step == 0)
    return;

  // This variable is used to compute the maximum number of particles
  // that can fit in the chosen insertion box before the acceptance function
  // is applied.
  unsigned int maximum_number_of_points = 1;

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
  // to the axis order.
  for (const unsigned int axis : axis_list)
    {
      AssertThrow(axis < dim,
                  ExcMessage("Insertion direction must be 0, 1 or 2"));

      axis_min[axis] = insertion_information.insertion_box_point_1(axis);
      axis_max[axis] = insertion_information.insertion_box_point_2(axis);

      // Assign the maximum number of insertion points according to the
      // direction and calculate the total number of points that fit in the box
      // (maximum_number_of_points = max_x * max_y * max_z)
      const unsigned int number_of_points =
        (axis_max[axis] - axis_min[axis]) /
        (insertion_information.distance_threshold * maximum_diameter);

      number_of_particles_directions[axis] = number_of_points;

      maximum_number_of_points *= number_of_points;
    }

  // We split the number of points evenly on every process.
  // The last process needs to be adjusted to match the maximum number of points
  // in the box.
  MPI_Comm communicator     = MPI_COMM_WORLD;
  auto     this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
  auto     n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);

  unsigned int n_points_this_proc =
    floor(maximum_number_of_points / n_mpi_process);
  if (this_mpi_process == (n_mpi_process - 1))
    n_points_this_proc =
      maximum_number_of_points -
      (n_mpi_process - 1) * floor(maximum_number_of_points / n_mpi_process);

  // First and last index of the unfiltered points that were assigned to this
  // process.
  unsigned int first_index_unfiltered;
  unsigned int last_index_unfiltered;
  if (this_mpi_process == (n_mpi_process - 1)) // For the last process
    {
      first_index_unfiltered = maximum_number_of_points - n_points_this_proc;
      last_index_unfiltered  = maximum_number_of_points;
    }
  else // For processes 1 to n-1
    {
      first_index_unfiltered = this_mpi_process * n_points_this_proc;
      last_index_unfiltered  = (this_mpi_process + 1) * n_points_this_proc;
    }

  // Now, we know that the indexes before applying the acceptance function
  // will go from first_index_unfiltered to last_index_unfiltered.
  // We count the number of insertion points that respect the acceptance
  // function and store the index of those points
  Point<dim> insertion_location;
  for (unsigned int index = first_index_unfiltered;
       index < last_index_unfiltered;
       ++index)
    {
      // Create the insertion point, associated with the current index, with
      // no offset
      find_insertion_location(
        insertion_location, index, 0., 0., insertion_information);

      // If the point respects the acceptance function, we insert the index in
      // the vector.
      if (acceptance_fct->value(insertion_location) > 0.)
        filted_box_index.push_back(index);
    }

  // Each process needs to know what is its first and last filtered indexes.
  // To do so, every process needs to know how many insertion point were kept
  // valid by the previous process.
  // The numeration used for the filtered indexes follows the same principal as
  // the one used before filtering. The only difference is that we skip the
  // rejected points.
  std::vector<unsigned int> starting_index_on_every_proc(n_mpi_process);
  starting_index_on_every_proc.resize(n_mpi_process);

  // The number of insertion points assigned to this process that respect the
  // acceptance function.
  const unsigned int n_valid_insertion_point_this_proc =
    filted_box_index.size();

  // Every process sends to process 0 the number of valid insertion points it
  // has.
  auto number_of_insertion_point_per_core =
    Utilities::MPI::gather(communicator, n_valid_insertion_point_this_proc, 0);

  if (this_mpi_process == 0)
    {
      starting_index_on_every_proc[0] = 0;
      for (unsigned int i = 1; i < n_mpi_process; ++i)
        starting_index_on_every_proc[i] =
          starting_index_on_every_proc[i - 1] +
          number_of_insertion_point_per_core[i - 1];

      number_of_valid_insertion_point_global =
        std::accumulate(number_of_insertion_point_per_core.begin(),
                        number_of_insertion_point_per_core.end(),
                        0);
    }

  // This scatters the number of particles to insert
  first_index_this_proc =
    Utilities::MPI::scatter(communicator, starting_index_on_every_proc, 0);

  // We find the last index with the size of the vector.
  last_index_this_proc =
    first_index_this_proc + n_valid_insertion_point_this_proc;

  // We need this for the adjust_insertion_count_by_insertion_box_capacity
  // function
  number_of_valid_insertion_point_global =
    Utilities::MPI::broadcast(communicator,
                              number_of_valid_insertion_point_global,
                              0);
}

template <int dim, typename PropertiesIndex>
void
InsertionVolume<dim, PropertiesIndex>::
  adjust_insertion_count_by_insertion_box_capacity(
    const InsertionInfo<dim> &insertion_information,
    const ConditionalOStream &pcout,
    unsigned int             &inserted_this_step)
{
  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step >
      number_of_valid_insertion_point_global)
    {
      pcout << "Warning: the requested number of particles for insertion ("
            << insertion_information.inserted_this_step
            << ") is higher than maximum expected number of particles ("
            << number_of_valid_insertion_point_global << ")" << std::endl;

      // Updating the number of inserted particles at each step
      inserted_this_step = number_of_valid_insertion_point_global;
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
