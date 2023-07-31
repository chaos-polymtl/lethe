#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/plane_insertion.h>

using namespace DEM;

// The constructor of plane insertion class.
template <int dim>
PlaneInsertion<dim>::PlaneInsertion(
  const DEMSolverParameters<dim> &                 dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation)
  : remained_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;
  particle_counter                = 0;

  // Initializing the cell that are close to the plan
  find_inplane_cells(
    triangulation,
    dem_parameters.insertion_info.insertion_plane_point,
    dem_parameters.insertion_info.insertion_plane_normal_vector);

  // Initializing the cell centers
  find_centers_of_inplane_cells();
}

template <int dim>
void
PlaneInsertion<dim>::find_inplane_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  Point<3>                                         plane_point,
  Tensor<1, 3>                                     plane_normal_vector)
{ // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by owned by the processor
      if (cell->is_locally_owned())
        {
          bool cell_in_plane = false;
          // Initializing the values

          // Vector that goes from the plane reference point to the vertex 0.
          Tensor<1, 3> connecting_vector_ref =
            point_nd_to_3d(cell->vertex(0)) - plane_point;

          // Normal distance between vertex 0 and the plane.
          double vertex_wall_distance_ref =
            connecting_vector_ref * plane_normal_vector;

          // Loop over all the vertices of the cell
          for (unsigned int vertex_id = 1; vertex_id < cell->n_vertices();
               ++vertex_id)
            {
              // Vector that goes from the plane reference point to the vertex
              // n-th
              Tensor<1, 3> connecting_vector =
                point_nd_to_3d(cell->vertex(vertex_id)) - plane_point;

              // Normal distance between vertex n-th and the plan
              double vertex_wall_distance =
                connecting_vector * plane_normal_vector;

              // If the multiplication of the two distances is negative, then at
              // least one of vertex of the cell is on both sides of the plane.
              // If the multiplication is equal to zero, then one vertex is on
              // the plane.
              if (vertex_wall_distance_ref * vertex_wall_distance <= 0)
                {
                  cell_in_plane = true;
                  break; // If so, we brake the loop, and go to the next cell
                }
            }
          if (cell_in_plane)
            {
              plane_cells_for_insertion.insert(cell);
            }
        }
    }
}

template <int dim>
void
PlaneInsertion<dim>::find_centers_of_inplane_cells()
{
  for (const auto &cell : plane_cells_for_insertion)
    {
      cells_centers.insert({cell->active_cell_index(), cell->center()});
    }
}

// The main insertion function. insert_particle is utilized to insert the
// particle at the cell center with a random shifts particles
template <int dim>
void
PlaneInsertion<dim>::insert(
  Particles::ParticleHandler<dim> &                particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim> &                 dem_parameters)
{
  if (remained_particles_of_each_type == 0 &&
      this->current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remained_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }
  if (remained_particles_of_each_type > 0)
    {
      MPI_Comm           communicator = triangulation.get_communicator();
      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);

      // List of the empty cells at this time step
      std::set<typename Triangulation<dim>::active_cell_iterator>
        empty_cells_on_proc;

      // Loop over the cell inplane
      for (const auto &cell : plane_cells_for_insertion)
        {
          // if the cell is empty ...
          if (particle_handler.n_particles_in_cell(cell) == 0)
            {
              // ... we add the cell in std::set
              empty_cells_on_proc.insert(cell);
            }
        }

      // Every processor knows how many particle (or empty cells) there is at
      // this time step. This is basically the only information share between
      // processor.
      auto number_of_particles_to_insert_per_core =
        Utilities::MPI::all_gather(communicator, empty_cells_on_proc.size());

      // Every processor does the following calculation. This could be done by a
      // "gather and scatter", but it's probably less efficient. First, we loop
      // over the vector that we have created using the all_gather. Every
      // integer in that vector gives the number of empty_cell of a certain
      // processor.
      for (unsigned int i = 0;
           i < number_of_particles_to_insert_per_core.size();
           ++i)
        {
          // We subtract those integer.
          remained_particles_of_each_type -=
            number_of_particles_to_insert_per_core[i];

          // If the remained_particles_of_each_type become negative, this means
          // that we are trying to insert to many particle in the simulation.
          // Thus, we need to decrease the number of empty_cell at this time
          // step. So...
          if (remained_particles_of_each_type <= 0)
            {
              // ... we decrease the number of particle at the current
              // processor, since it's the one who got the variable to become
              // negative. Then...
              number_of_particles_to_insert_per_core[i] +=
                remained_particles_of_each_type;

              // ... we put remained_particles_of_each_type back to zero.
              remained_particles_of_each_type = 0;

              // Loop over the remaining processor and put those at zero.
              for (unsigned int j = i + 1;
                   j < number_of_particles_to_insert_per_core.size();
                   ++j)
                {
                  number_of_particles_to_insert_per_core[j] = 0;
                }
              break;
            }
        }

      // Now, every processor knows how many particle it needs to insert. They
      // also know how many particle the other processor have to insert. This
      // way, it's possible for a processor to find the correct range of
      // particle id to use with the particle_counter variable.

      // This is the starting id of the processor.
      unsigned int id_on_proc =
        std::accumulate(number_of_particles_to_insert_per_core.begin(),
                        number_of_particles_to_insert_per_core.begin() +
                          this_mpi_process,
                        particle_counter);

      // We yet didn't say which empty cells won't be use anymore.
      // empty_cells_on_proc needs to be shorter for some processor.
      while (empty_cells_on_proc.size() <
             number_of_particles_to_insert_per_core[this_mpi_process])
        {
          auto it = empty_cells_on_proc.end();
          empty_cells_on_proc.erase(it);
        }

      // There's probably a better ways to define the properties of a particle,
      // but for now it works...
      double type     = this->current_inserting_particle_type;
      double diameter = dem_parameters.lagrangian_physical_properties
                          .particle_average_diameter.at(type);
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

      // Loop over the empty cells we kept.
      for (const auto &cell : empty_cells_on_proc)
        {
          // Position of the center of the cell
          Point<dim> insertion_location =
            cells_centers.at(cell->active_cell_index());

          // Generate the random point of insertion
          insertion_location(0) +=
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX) *
            dem_parameters.insertion_info.random_number_range;

          insertion_location(1) +=
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX) *
            dem_parameters.insertion_info.random_number_range;

          if constexpr (dim == 3)
            {
              insertion_location(2) +=
                static_cast<double>(rand()) / static_cast<double>(RAND_MAX) *
                dem_parameters.insertion_info.random_number_range;
            }

          // Insertion
          Point<dim> ref_point;
          particle_handler.insert_particle(insertion_location,
                                           ref_point,
                                           id_on_proc++,
                                           cell,
                                           properties_of_one_particle);
        }
      // We need to synchronize the particle_counter on all the processor for
      // the next insertion time step.
      particle_counter +=
        std::accumulate(number_of_particles_to_insert_per_core.begin(),
                        number_of_particles_to_insert_per_core.end(),
                        0);
    }
}
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
