#include <dem/non_uniform_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (nx, ny and nz) are also obtained.
template <int dim>
NonUniformInsertion<dim>::NonUniformInsertion(
  const DEMSolverParameters<dim> &dem_parameters,
  unsigned int &                  inserted_this_step,
  unsigned int &                  nx,
  unsigned int &                  ny,
  unsigned int &                  nz)
{
  // Getting properties as local parameters
  const auto physical_properties   = dem_parameters.physicalProperties;
  const auto insertion_information = dem_parameters.insertionInfo;

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

  // nx, ny and nz are the results of discretization of the insertion domain in
  // x, y and z directions
  nx = int(
    (insertion_information.x_max - insertion_information.x_min) /
    (insertion_information.distance_threshold * physical_properties.diameter));
  ny = int(
    (insertion_information.y_max - insertion_information.y_min) /
    (insertion_information.distance_threshold * physical_properties.diameter));
  if (dim == 3)
    {
      nz = int((insertion_information.z_max - insertion_information.z_min) /
               (insertion_information.distance_threshold *
                physical_properties.diameter));
    }
}

// The main insertion function. Insert_global_function is utilized to insert the
// particles
template <int dim>
void
NonUniformInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters,
  unsigned int &                                   inserted_this_step,
  const unsigned int &                             nx,
  const unsigned int &                             ny,
  const unsigned int &                             nz,
  unsigned int &                                   remained_particles)
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
      insertion_points = this->assign_insertion_points(
        dem_parameters, inserted_this_step, nx, ny, nz);

      // Assigning inserted particles properties using
      // assign_particle_properties function
      std::vector<std::vector<double>> particle_properties;
      particle_properties =
        this->assign_particle_propertis(dem_parameters,
                                        inserted_this_step,
                                        particle_handler.n_global_particles());

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points,
                                               global_bounding_boxes,
                                               particle_properties);

      // Updating remaining particles
      remained_particles -= inserted_this_step;
    }
}

// This function creates a vector of random doubles using the input paramteres
// in the paramter handler
template <int dim>
std::vector<double>
NonUniformInsertion<dim>::create_random_number_container(
  const int &   inserted_this_step,
  const double &random_number_range,
  const int &   random_number_seed)
{
  std::vector<double> random_container;
  for (int i = 0; i < inserted_this_step; ++i)
    {
      srand(random_number_seed * (i + 1));
      random_container.push_back((((double)rand()) / ((double)RAND_MAX)) *
                                 random_number_range);
    }
  return random_container;
}

// This function assigns the insertion points of the inserted particles
template <int dim>
std::vector<Point<dim>>
NonUniformInsertion<dim>::assign_insertion_points(
  const DEMSolverParameters<dim> &dem_parameters,
  const unsigned int &            inserted_this_step,
  const unsigned int &            nx,
  const unsigned int &            ny,
  const unsigned int &            nz)
{
  // Initilizing the output vector
  std::vector<Point<dim>> insertion_positions;

  // Getting properties as local parameters
  const auto physical_properties   = dem_parameters.physicalProperties;
  const auto insertion_information = dem_parameters.insertionInfo;

  // Calling random number generator
  std::vector<double> random_number_vector;
  random_number_vector = this->create_random_number_container(
    inserted_this_step,
    insertion_information.random_number_range,
    insertion_information.random_number_seed);

  // Creating a particle counter
  unsigned int particle_counter = 0;

  for (unsigned int i = 0; i < nx; ++i)
    for (unsigned int j = 0; j < ny; ++j)
      {
        // Adapt the last index to the dimensionality of the problem
        unsigned int dim_nz = (dim == 3) ? nz : 1;
        for (unsigned int k = 0; k < dim_nz; ++k)
          {
            // We need to check if the number of inserted particles so far at
            // this step (particle_counter) reached the total desired number of
            // inserted particles at this step
            if (particle_counter < inserted_this_step)
              {
                Point<dim> position;
                // Obtaning position of the inserted particle
                // In non-uniform insertion, random numbers were created and are
                // added to the position of particles. In order to create more
                // randomness, the random vector is read once from the beginning
                // and once from the end to be used in positions [0] and [1]
                position[0] = insertion_information.x_min +
                              (physical_properties.diameter / 2) +
                              (i * insertion_information.distance_threshold *
                               physical_properties.diameter) +
                              random_number_vector[particle_counter] *
                                physical_properties.diameter;
                position[1] = insertion_information.y_min +
                              (physical_properties.diameter / 2) +
                              (j * insertion_information.distance_threshold *
                               physical_properties.diameter) +
                              random_number_vector[insertion_information
                                                     .inserted_this_step -
                                                   particle_counter] *
                                physical_properties.diameter;
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

template class NonUniformInsertion<2>;
template class NonUniformInsertion<3>;
