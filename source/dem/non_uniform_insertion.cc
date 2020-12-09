#include <deal.II/base/point.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <dem/dem_properties.h>
#include <dem/non_uniform_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim>
NonUniformInsertion<dim>::NonUniformInsertion(
  const DEMSolverParameters<dim> &dem_parameters,
  const double &                  maximum_particle_diameter)
  : remained_particles_of_each_type(
      dem_parameters.physical_properties.number.at(0))
{
  // Inializing current inserting particle type
  current_inserting_particle_type = 0;

  this->inserted_this_step = 0;

  this->maximum_diameter = maximum_particle_diameter;
}

// The main insertion function. Insert_global_function is utilized to insert the
// particles
template <int dim>
void
NonUniformInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters)
{
  if (remained_particles_of_each_type == 0 &&
      current_inserting_particle_type !=
        dem_parameters.physical_properties.particle_type_number - 1)
    {
      remained_particles_of_each_type =
        dem_parameters.physical_properties.number.at(
          ++current_inserting_particle_type);
    }

  // Check to see if the remained uninserted particles is equal to zero or not
  if (remained_particles_of_each_type != 0)
    {
      this->calculate_insertion_domain_maximum_particle_number(dem_parameters);

      // The inserted_this_step value is the mimnum of
      // remained_particles_of_each_type and inserted_this_step
      this->inserted_this_step =
        std::min(remained_particles_of_each_type, this->inserted_this_step);

      MPI_Comm communicator = triangulation.get_communicator();

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      unsigned int this_mpi_process =
        Utilities::MPI::this_mpi_process(communicator);

      // Finding insertion points using assign_insertion_points function
      std::vector<Point<dim>> insertion_points;
      insertion_points.resize(0);
      if (this_mpi_process == 0)
        insertion_points =
          this->assign_insertion_points(dem_parameters.insertion_info);

      // Assigning inserted particles properties using
      // assign_particle_properties function
      if (this_mpi_process == 0)
        this->assign_particle_properties(dem_parameters,
                                         this->inserted_this_step,
                                         current_inserting_particle_type,
                                         this->particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points,
                                               global_bounding_boxes,
                                               this->particle_properties);

      // Updating remaining particles
      remained_particles_of_each_type -= this->inserted_this_step;

      this->print_insertion_info(this->inserted_this_step,
                                 remained_particles_of_each_type,
                                 current_inserting_particle_type);
    }
}

// This function creates a vector of random doubles using the input paramteres
// in the paramter handler
template <int dim>
std::vector<double>
NonUniformInsertion<dim>::create_random_number_container(
  const double &random_number_range,
  const int &   random_number_seed)
{
  std::vector<double> random_container;
  for (unsigned int i = 0; i < this->inserted_this_step; ++i)
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
  const Parameters::Lagrangian::InsertionInfo &insertion_information)
{
  // Initilizing the output vector
  std::vector<Point<dim>> insertion_positions;

  // Calling random number generator
  std::vector<double> random_number_vector;
  random_number_vector = this->create_random_number_container(
    insertion_information.random_number_range,
    insertion_information.random_number_seed);

  // Creating a particle counter
  unsigned int particle_counter = 0;

  for (unsigned int i = 0; i < this->number_of_particles_x_direction; ++i)
    for (unsigned int j = 0; j < this->number_of_particles_y_direction; ++j)
      {
        // Adapt the last index to the dimensionality of the problem
        unsigned int dim_nz =
          (dim == 3) ? this->number_of_particles_z_direction : 1;
        for (unsigned int k = 0; k < dim_nz; ++k)
          {
            // We need to check if the number of inserted particles so far at
            // this step (particle_counter) reached the total desired number of
            // inserted particles at this step
            if (particle_counter < this->inserted_this_step)
              {
                Point<dim> position;
                // Obtaning position of the inserted particle
                // In non-uniform insertion, random numbers were created and are
                // added to the position of particles. In order to create more
                // randomness, the random vector is read once from the beginning
                // and once from the end to be used in positions [0] and [1]
                position[0] = insertion_information.x_min +
                              (this->maximum_diameter / 2) +
                              (i * insertion_information.distance_threshold *
                               this->maximum_diameter) +
                              random_number_vector[particle_counter] *
                                this->maximum_diameter;
                position[1] = insertion_information.y_min +
                              (this->maximum_diameter / 2) +
                              (j * insertion_information.distance_threshold *
                               this->maximum_diameter) +
                              random_number_vector[this->inserted_this_step -
                                                   particle_counter - 1] *
                                this->maximum_diameter;
                if (dim == 3)
                  {
                    position[2] =
                      insertion_information.z_min +
                      (this->maximum_diameter / 2) +
                      (k * insertion_information.distance_threshold *
                       this->maximum_diameter);
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
