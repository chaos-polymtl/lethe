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

  this->calculate_insertion_domain_maximum_particle_number(dem_parameters);
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
      this->calculate_insertion_domain_maximum_particle_number(dem_parameters);

      // The inserted_this_step value is the mimnum of remained_particles and
      // inserted_this_step
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

// This function assigns the insertion points of the inserted particles
template <int dim>
std::vector<Point<dim>>
UniformInsertion<dim>::assign_insertion_points(
  const Parameters::Lagrangian::InsertionInfo &insertion_information)
{
  std::vector<Point<dim>> insertion_positions;

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
            // Check if the number of inserted particles so far at
            // this step reached the total desired number of inserted particles
            // at this step
            if (particle_counter < this->inserted_this_step)
              {
                Point<dim> position;
                // Obtaning position of the inserted particle
                position[0] = insertion_information.x_min +
                              (this->maximum_diameter / 2) +
                              (i * insertion_information.distance_threshold *
                               this->maximum_diameter);
                position[1] = insertion_information.y_min +
                              (this->maximum_diameter / 2) +
                              (j * insertion_information.distance_threshold *
                               this->maximum_diameter);
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

template class UniformInsertion<2>;
template class UniformInsertion<3>;
