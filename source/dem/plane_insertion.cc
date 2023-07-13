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
  this->insertion_plane_normal_vector = dem_parameters.insertion_info.insertion_plane_normal_vector;
  this->insertion_plane_point = dem_parameters.insertion_info.insertion_plane_point;
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

template <>
void
PlaneInsertion<3>::find_plane_cells_for_plane_insertion(
  const parallel::distributed::Triangulation<3> &triangulation)
{
  const double maximal_cell_diameter =
    GridTools::maximal_cell_diameter(triangulation);

  // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by owned by the processor
      if (cell->is_locally_owned())
        {
          bool cell_in_plane = true;
          for (unsigned int vertex_id = 0; vertex_id < cell->n_vertices(); ++vertex_id)
          {
            Tensor<1,3> connecting_vector =
            cell->vertex(vertex_id) - insertion_plane_point;
            double vertex_wall_distance = std::abs(connecting_vector * insertion_plane_normal_vector);

            if ( vertex_wall_distance > maximal_cell_diameter)
              {
                cell_in_plane = false;
                break;
              }
          }
          if (cell_in_plane)
          {
            plane_cells_for_insertion.insert(cell);
          }
        }
    }
}

template <>
void
PlaneInsertion<3>::find_center_of_in_plane_cells()
{
  for (const auto &cell : plane_cells_for_insertion)
    {

    }

}



;
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
