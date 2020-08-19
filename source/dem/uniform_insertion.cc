#include <dem/uniform_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim>
UniformInsertion<dim>::UniformInsertion(
  const DEMSolverParameters<dim> &dem_parameters)
  : remained_particles(dem_parameters.insertion_info.total_particle_number)
  , inserted_this_step(0)
{
  // Getting properties as local parameters
  const auto physical_properties   = dem_parameters.physical_properties;
  const auto insertion_information = dem_parameters.insertion_info;

  // This variable is used for calculation of the maximum number of particles
  // that can fit in the chosen insertion box
  int maximum_particle_number;

  // distance_threshold shows the ratio of the distance between the centers of
  // two adjacent particles to the diameter of particles
  maximum_particle_number =
    int((insertion_information.x_max - insertion_information.x_min) /
        (insertion_information.distance_threshold *
         physical_properties.diameter)) *
    int((insertion_information.y_max - insertion_information.y_min) /
        (insertion_information.distance_threshold *
         physical_properties.diameter));
  if (dim == 3)
    {
      maximum_particle_number =
        maximum_particle_number *
        int((insertion_information.z_max - insertion_information.z_min) /
            (insertion_information.distance_threshold *
             physical_properties.diameter));
    }

  // If the inserted number of particles at this step exceeds the maximum
  // number, a warning is printed
  if (insertion_information.inserted_this_step > maximum_particle_number)
    {
      std::cout << "The inserted number of particles ("
                << insertion_information.inserted_this_step
                << ") is higher than maximum expected number of particles ("
                << maximum_particle_number << ")" << std::endl;
      std::cout << "Inserting " << maximum_particle_number << " at this step"
                << std::endl;

      // Updating the number of inserted particles at each step
      inserted_this_step = maximum_particle_number;
    }
  else
    {
      inserted_this_step = insertion_information.inserted_this_step;
    }

  // number_of_particles_x_direction, number_of_particles_y_direction and
  // number_of_particles_z_direction are the results of discretization of the
  // insertion domain in x, y and z directions
  number_of_particles_x_direction = int(
    (insertion_information.x_max - insertion_information.x_min) /
    (insertion_information.distance_threshold * physical_properties.diameter));
  number_of_particles_y_direction = int(
    (insertion_information.y_max - insertion_information.y_min) /
    (insertion_information.distance_threshold * physical_properties.diameter));
  if (dim == 3)
    {
      number_of_particles_z_direction =
        int((insertion_information.z_max - insertion_information.z_min) /
            (insertion_information.distance_threshold *
             physical_properties.diameter));
    }
}

// The main insertion function. Insert_global_function is utilized to insert the
// particles
template <int dim>
void
UniformInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters)
{
  // Check to see if the remained uninserted particles is equal to zero or not
  if (remained_particles != 0)
    {
      // The inserted_this_step value is the mimnum of remained_particles and
      // inserted_this_step
      inserted_this_step = std::min(remained_particles, inserted_this_step);

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

      // Finding insertion points using assign_insertion_points function
      std::vector<Point<dim>> insertion_points;
      insertion_points = this->assign_insertion_points(dem_parameters);

      // Assigning inserted particles properties using
      // assign_particle_properties function
      std::vector<std::vector<double>> particle_properties;
      particle_properties =
        this->assign_particle_properties(dem_parameters,
                                         inserted_this_step,
                                         particle_handler.n_global_particles());

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points,
                                               global_bounding_boxes,
                                               particle_properties);

      // Updating remaining particles
      remained_particles -= inserted_this_step;

      this->print_insertion_info(inserted_this_step, remained_particles);
    }
}

// This function assigns the insertion points of the inserted particles
template <int dim>
std::vector<Point<dim>>
UniformInsertion<dim>::assign_insertion_points(
  const DEMSolverParameters<dim> &dem_parameters)
{
  std::vector<Point<dim>> insertion_positions;

  // Getting properties as local parameters
  const auto physical_properties   = dem_parameters.physical_properties;
  const auto insertion_information = dem_parameters.insertion_info;

  // Creating a particle counter
  unsigned int particle_counter = 0;

  for (unsigned int i = 0; i < number_of_particles_x_direction; ++i)
    for (unsigned int j = 0; j < number_of_particles_y_direction; ++j)
      {
        // Adapt the last index to the dimensionality of the problem
        unsigned int dim_nz = (dim == 3) ? number_of_particles_z_direction : 1;
        for (unsigned int k = 0; k < dim_nz; ++k)
          {
            // We need to check if the number of inserted particles so far at
            // this step reached the total desired number of inserted particles
            // at this step
            if (particle_counter < inserted_this_step)
              {
                Point<dim> position;
                // Obtaning position of the inserted particle
                position[0] = insertion_information.x_min +
                              (physical_properties.diameter / 2) +
                              (i * insertion_information.distance_threshold *
                               physical_properties.diameter);
                position[1] = insertion_information.y_min +
                              (physical_properties.diameter / 2) +
                              (j * insertion_information.distance_threshold *
                               physical_properties.diameter);
                // Adding a threshold distance to even rows of insertion
                if (k % 2 == 0)
                  {
                    position[0] =
                      position[0] + (physical_properties.diameter) / 2.0;
                    position[1] =
                      position[1] + (physical_properties.diameter) / 2.0;
                  }
                if (dim == 3)
                  {
                    position[2] =
                      insertion_information.z_min +
                      (physical_properties.diameter / 2) +
                      (k * insertion_information.distance_threshold *
                       physical_properties.diameter);
                  }

                insertion_positions.push_back(position);
                particle_counter++;
              }
          }
      }
  return insertion_positions;
}

template class UniformInsertion<2>;
template class UniformInsertion<3>;
