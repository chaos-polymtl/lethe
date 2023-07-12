#include <dem/plane_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim>
PlaneInsertion<dim>::PlaneInsertion(
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
PlaneInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters)
{

}

// This function creates a vector of random doubles using the input paramteres
// in the parameter handler
template <int dim>
void
PlaneInsertion<dim>::create_random_number_container(
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
void
PlaneInsertion<3>::find_insertion_location_plane(
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
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
