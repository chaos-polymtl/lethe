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

      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      auto n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);

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

      // Distbuting particles between processors
      if (this_mpi_process != (n_mpi_process - 1))
        this->inserted_this_step_this_proc =
          floor(this->inserted_this_step / n_mpi_process);
      else
        this->inserted_this_step_this_proc =
          this->inserted_this_step -
          (n_mpi_process - 1) * floor(this->inserted_this_step / n_mpi_process);

      // Finding insertion points using assign_insertion_points function
      std::vector<Point<dim>> global_insertion_points;
      std::vector<Point<dim>> insertion_points_on_proc;
      global_insertion_points.reserve(this->inserted_this_step);
      insertion_points_on_proc.reserve(this->inserted_this_step_this_proc);

      this->assign_insertion_points(global_insertion_points,
                                    dem_parameters.insertion_info);

      // Distributing the insertion points between the processors
      if (this_mpi_process != (n_mpi_process - 1))
        {
          int begin_element =
            this_mpi_process * this->inserted_this_step_this_proc;
          int end_element =
            (this_mpi_process + 1) * this->inserted_this_step_this_proc - 1;
          insertion_points_on_proc = std::vector<Point<dim>>(
            global_insertion_points.cbegin() + begin_element,
            global_insertion_points.cbegin() + end_element + 1);
        }
      else
        {
          int begin_element =
            this->inserted_this_step - this->inserted_this_step_this_proc;
          insertion_points_on_proc =
            std::vector<Point<dim>>(global_insertion_points.cbegin() +
                                      begin_element,
                                    global_insertion_points.cend());
        }

      // Assigning inserted particles properties using
      // assign_particle_properties function
      this->assign_particle_properties(dem_parameters,
                                       this->inserted_this_step_this_proc,
                                       current_inserting_particle_type,
                                       this->particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(insertion_points_on_proc,
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
  const Parameters::Lagrangian::InsertionInfo &insertion_information)
{
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
            if (particle_counter < this->inserted_this_step)
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
                    position[0] = position[0] + (this->maximum_diameter) / 2.0;
                    position[1] = position[1] + (this->maximum_diameter) / 2.0;
                  }
                if (dim == 3)
                  {
                    position[2] =
                      insertion_information.z_min +
                      +((k + 0.5) * insertion_information.distance_threshold) *
                        this->maximum_diameter;
                  }

                insertion_positions.push_back(position);
                particle_counter++;
              }
          }
      }
}

template class UniformInsertion<2>;
template class UniformInsertion<3>;
