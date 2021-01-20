#include <dem/uniform_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim>
UniformInsertion<dim>::UniformInsertion(
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

// The main insertion function. Insert_global_function is used to insert the
// particles
template <int dim>
void
UniformInsertion<dim>::insert(
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
      MPI_Comm           communicator = triangulation.get_communicator();
      ConditionalOStream pcout        = {
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0};

      this->calculate_insertion_domain_maximum_particle_number(dem_parameters,
                                                               pcout);

      // The inserted_this_step value is the mimnum of remained_particles and
      // inserted_this_step
      this->inserted_this_step =
        std::min(remained_particles_of_each_type, this->inserted_this_step);

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // Finding insertion points using assign_insertion_points function
      std::vector<Point<dim>> insertion_points;
      insertion_points.reserve(this->inserted_this_step);
      this->assign_insertion_points(insertion_points,
                                    dem_parameters.insertion_info,
                                    communicator);

      // Assigning inserted particles properties using
      // assign_particle_properties function
      this->particle_properties.reserve(this->inserted_this_step);
      this->assign_particle_properties(dem_parameters,
                                       this->inserted_this_step,
                                       current_inserting_particle_type,
                                       this->particle_properties,
                                       communicator);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points,
                                               global_bounding_boxes,
                                               this->particle_properties);

      // Updating remaining particles
      remained_particles_of_each_type -= this->inserted_this_step;

      this->print_insertion_info(this->inserted_this_step,
                                 remained_particles_of_each_type,
                                 current_inserting_particle_type,
                                 pcout);
    }
}

// This function assigns the insertion points of the inserted particles
template <int dim>
void
UniformInsertion<dim>::assign_insertion_points(
  std::vector<Point<dim>> &                    insertion_positions,
  const Parameters::Lagrangian::InsertionInfo &insertion_information,
  const MPI_Comm &                             communicator)
{
  unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(communicator);

  // Creating a particle counter
  unsigned int particle_counter     = 0;
  unsigned int particle_counter_sum = 0;

  for (unsigned int i = 0; i < this->number_of_particles_x_direction; ++i)
    for (unsigned int j = 0; j < this->number_of_particles_y_direction; ++j)
      {
        // Adapt the last index to the dimensionality of the problem
        unsigned int dim_nz =
          (dim == 3) ? this->number_of_particles_z_direction : 1;
        for (unsigned int k = 0; k < dim_nz; ++k)
          {
            if (this->particle_on_processor(i,
                                            j,
                                            k,
                                            this_mpi_process,
                                            Utilities::MPI::n_mpi_processes(
                                              communicator)))
              {
                // Check if the number of inserted particles so far at
                // this step reached the total desired number of inserted
                // particles at this step
                particle_counter_sum =
                  Utilities::MPI::sum(particle_counter_sum, communicator);
                if (particle_counter_sum < this->inserted_this_step)
                  {
                    Point<dim> position;
                    // Obtaning position of the inserted particle
                    position[0] =
                      insertion_information.x_min +
                      ((i + 0.5) * insertion_information.distance_threshold) *
                        this->maximum_diameter;
                    position[1] =
                      insertion_information.y_min +
                      +((j + 0.5) * insertion_information.distance_threshold) *
                        this->maximum_diameter;

                    // Adding a threshold distance to even rows of insertion
                    if (k % 2 == 0)
                      {
                        position[0] =
                          position[0] + (this->maximum_diameter) / 2.0;
                        position[1] =
                          position[1] + (this->maximum_diameter) / 2.0;
                      }
                    if (dim == 3)
                      {
                        position[2] =
                          insertion_information.z_min +
                          +((k + 0.5) *
                            insertion_information.distance_threshold) *
                            this->maximum_diameter;
                      }

                    insertion_positions.push_back(position);
                    particle_counter++;
                  }
              }
          }
      }
}

template class UniformInsertion<2>;
template class UniformInsertion<3>;
