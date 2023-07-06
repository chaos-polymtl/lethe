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
  const double                    maximum_particle_diameter)
  : remained_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
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
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remained_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }

  // Check to see if the remained uninserted particles is equal to zero or not
  if (remained_particles_of_each_type != 0)
    {
      MPI_Comm           communicator = triangulation.get_communicator();
      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      auto n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);


      this->calculate_insertion_domain_maximum_particle_number(dem_parameters,
                                                               pcout);
      // The inserted_this_step value is the mimnum of
      // remained_particles_of_each_type and inserted_this_step
      this->inserted_this_step =
        std::min(remained_particles_of_each_type, this->inserted_this_step);

      // Obtaining global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // Distbuting particles between processors
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
        dem_parameters.insertion_info.random_number_range,
        dem_parameters.insertion_info.random_number_seed);

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
          find_insertion_location_nonuniform(
            insertion_location,
            id,
            random_number_vector[particle_counter],
            random_number_vector[this->inserted_this_step - particle_counter -
                                 1],
            dem_parameters.insertion_info);
          insertion_points_on_proc.push_back(insertion_location);
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

// This function creates a vector of random doubles using the input paramteres
// in the parameter handler
template <int dim>
void
NonUniformInsertion<dim>::create_random_number_container(
  std::vector<double> &random_container,
  const double         random_number_range,
  const int            random_number_seed)
{
  for (unsigned int i = 0; i < this->inserted_this_step; ++i)
    {
      srand(random_number_seed * (i + 1));
      random_container.push_back((((double)rand()) / ((double)RAND_MAX)) *
                                 random_number_range);
    }
}

// This function assigns the insertion points of the inserted particles
template <>
void NonUniformInsertion<2>::find_insertion_location_nonuniform(
  Point<2> &                                   insertion_location,
  const unsigned int                           id,
  const double                                 random_number1,
  const double                                 random_number2,
  const Parameters::Lagrangian::InsertionInfo &insertion_information)
{
  std::vector<int> insertion_index;
  insertion_index.resize(2);

  unsigned int axis_0, axis_1;
  int          number_of_particles_0;

  // First direction (axis) to have particles inserted
  axis_0                  = insertion_information.axis_0;
  number_of_particles_0   = this->number_of_particles_directions[axis_0];
  insertion_index[axis_0] = id % number_of_particles_0;
  insertion_location[axis_0] =
    this->axis_min[axis_0] + ((insertion_index[axis_0] + 0.5) *
                                insertion_information.distance_threshold -
                              random_number1) *
                               this->maximum_diameter;

  // Second direction (axis) to have particles inserted
  axis_1                  = insertion_information.axis_1;
  insertion_index[axis_1] = static_cast<int>(id / number_of_particles_0);
  insertion_location[axis_1] =
    this->axis_min[axis_1] + ((insertion_index[axis_1] + 0.5) *
                                insertion_information.distance_threshold -
                              random_number2) *
                               this->maximum_diameter;
}

template <>
void NonUniformInsertion<3>::find_insertion_location_nonuniform(
  Point<3> &                                   insertion_location,
  const unsigned int                           id,
  const double                                 random_number1,
  const double                                 random_number2,
  const Parameters::Lagrangian::InsertionInfo &insertion_information)
{
  std::vector<int> insertion_index;
  insertion_index.resize(3);

  unsigned int axis_0, axis_1, axis_2;
  int          number_of_particles_0, number_of_particles_1;

  // First direction (axis) to have particles inserted
  axis_0                  = insertion_information.axis_0;
  number_of_particles_0   = this->number_of_particles_directions[axis_0];
  insertion_index[axis_0] = id % number_of_particles_0;
  insertion_location[axis_0] =
    this->axis_min[axis_0] + ((insertion_index[axis_0] + 0.5) *
                                insertion_information.distance_threshold -
                              random_number1) *
                               this->maximum_diameter;

  // Second direction (axis) to have particles inserted
  axis_1                = insertion_information.axis_1;
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
  axis_2 = insertion_information.axis_2;
  insertion_index[axis_2] =
    static_cast<int>(id / (number_of_particles_0 * number_of_particles_1));
  insertion_location[axis_2] =
    this->axis_min[axis_2] + ((insertion_index[axis_2] + 0.5) *
                                insertion_information.distance_threshold -
                              random_number1) *
                               this->maximum_diameter;
}

template class NonUniformInsertion<2>;
template class NonUniformInsertion<3>;
