#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/plane_insertion.h>

using namespace DEM;

// The constructor of non-uniform insertion class. In the constructor, we
// investigate if the insertion box is adequately large to handle the desired
// number of inserted particles. The number of insertion points in each
// direction (number_of_particles_x_direction, number_of_particles_y_direction
// and number_of_particles_z_direction) are also obtained
template <int dim>
PlaneInsertion<dim>::PlaneInsertion(
  const DEMSolverParameters<dim> &                 dem_parameters,
  const double                                     maximum_particle_diameter,
  const parallel::distributed::Triangulation<dim> &triangulation)
{


  this->inserted_this_step = 0;
  this->maximum_diameter   = maximum_particle_diameter;
  this->insertion_plane_normal_vector =
    dem_parameters.insertion_info.insertion_plane_normal_vector;
  this->insertion_plane_point =
    dem_parameters.insertion_info.insertion_plane_point;

  current_inserting_particle_type = 0;
  particle_counter                = 0;
  minimal_cell_diameter = GridTools::minimal_cell_diameter(triangulation) / std::sqrt(3.);
  srand(dem_parameters.insertion_info.random_number_seed);

  find_plane_cells_for_plane_insertion(triangulation);
  find_center_of_inplane_cells();
}

// This function creates a vector of random doubles using the input parameters
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

template <int dim>
void
PlaneInsertion<dim>::find_plane_cells_for_plane_insertion(
  const parallel::distributed::Triangulation<dim> &triangulation)
{ // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by owned by the processor
      if (cell->is_locally_owned())
        {
          bool cell_in_plane = true;
          // Loop over all the vertices of the cell
          for (unsigned int vertex_id = 0; vertex_id < cell->n_vertices();
               ++vertex_id)
            {
              // Vector that goes from the plane reference point to the vertice
              Tensor<1, 3> connecting_vector =
                point_nd_to_3d(cell->vertex(vertex_id)) - insertion_plane_point;

              // Normal distance of that vector
              double vertex_wall_distance =
                std::abs(connecting_vector * insertion_plane_normal_vector);

              // If the distance of one of the vertices if more than the
              // minima_cell_diameter, than the cell is not in the plane.
              if (vertex_wall_distance > (minimal_cell_diameter))
                {
                  cell_in_plane = false;
                  break; // If so, we brake the loop, and go to the next cell
                }
            }
          if (cell_in_plane) // Not sure if this if is necessary since we use a
                             // break
            {
              plane_cells_for_insertion.insert(cell);
            }
        }
    }
}

template <int dim>
void
PlaneInsertion<dim>::find_center_of_inplane_cells()
{
  for (const auto &cell : plane_cells_for_insertion)
    {
      cells_centers.insert({cell->active_cell_index(), cell->center()});
    }
}

// The main insertion function. insert_particle is utilized to insert the
// particle at the cell center with a rawdom shifts particles
template <int dim>
void
PlaneInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &particle_handler,
  const parallel::distributed::Triangulation<dim> & /*triangulation*/,
  const DEMSolverParameters<dim> &dem_parameters)
{
  double type = 0;
  double diameter =
    dem_parameters.lagrangian_physical_properties.particle_average_diameter.at(
      type);
  double density =
    dem_parameters.lagrangian_physical_properties.density_particle.at(type);
  double vel_x        = dem_parameters.insertion_info.vel_x;
  double vel_y        = dem_parameters.insertion_info.vel_y;
  double vel_z        = dem_parameters.insertion_info.vel_z;
  double omega_x      = dem_parameters.insertion_info.omega_x;
  double omega_y      = dem_parameters.insertion_info.omega_y;
  double omega_z      = dem_parameters.insertion_info.omega_z;
  double fem_force_x  = 0.;
  double fem_force_y  = 0.;
  double fem_force_z  = 0.;
  double fem_torque_x = 0.;
  double fem_torque_y = 0.;
  double fem_torque_z = 0.;
  double mass         = density * 4. / 3. * M_PI *
                Utilities::fixed_power<3, double>(diameter * 0.5);
  double volumetric_contribution = 0.;

  std::vector<double> properties_of_one_particle{type,
                                                 diameter,
                                                 vel_x,
                                                 vel_y,
                                                 vel_z,
                                                 omega_x,
                                                 omega_y,
                                                 omega_z,
                                                 fem_force_x,
                                                 fem_force_y,
                                                 fem_force_z,
                                                 fem_torque_x,
                                                 fem_torque_y,
                                                 fem_torque_z,
                                                 mass,
                                                 volumetric_contribution};
  // Loop over the cell inplane
  for (const auto &cell : plane_cells_for_insertion)
    {
     Point<dim> particle_location= cells_centers.at(
        cell->active_cell_index());
      Point<dim> ref_point;
      // if the cell is empty
      if (particle_handler.n_particles_in_cell(cell) == 0)
        {
                   particle_location(0) += static_cast<double>(rand()) /
                         static_cast<double>(RAND_MAX) *
                         dem_parameters.insertion_info.random_number_range;

          particle_location(1) += static_cast<double>(rand()) /
                         static_cast<double>(RAND_MAX) *
                         dem_parameters.insertion_info.random_number_range;

          if constexpr (dim == 3)
            {
              particle_location(2) += static_cast<double>(rand()) /
                             static_cast<double>(RAND_MAX) *
                             dem_parameters.insertion_info.random_number_range;
            }
          // insert a particle in the middle of the cell
          particle_handler.insert_particle(particle_location,
                                           ref_point,
                                           particle_counter++,
                                           cell,
                                           properties_of_one_particle);
        }
    }
}



;
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
