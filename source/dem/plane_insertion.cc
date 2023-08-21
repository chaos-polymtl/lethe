#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/plane_insertion.h>

using namespace DEM;

// The constructor of plane insertion class. In the constructor, we find which
// cells are going to be use for the insertion and we also find the centers of
// those cells.
template <int dim>
PlaneInsertion<dim>::PlaneInsertion(
  const DEMSolverParameters<dim> &                 dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation)
  : remained_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;

  // Finding which cells are inplane
  this->find_inplane_cells(
    triangulation,
    dem_parameters.insertion_info.insertion_plane_point,
    dem_parameters.insertion_info.insertion_plane_normal_vector);
  // // Finding the center of those cells
  this->find_centers_of_inplane_cells();
  mark_for_update = false;

  change_to_triangulation =
    triangulation.signals.any_change.connect([&] { mark_for_update = true; });
}

template <int dim>
void
PlaneInsertion<dim>::find_inplane_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  Point<3>                                         plane_point,
  Tensor<1, 3>                                     plane_normal_vector)
{ // Looping through cells
  plane_cells_for_insertion.clear();
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
  cells_centers.clear();
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
      if (mark_for_update)
        {
          find_inplane_cells(
            triangulation,
            dem_parameters.insertion_info.insertion_plane_point,
            dem_parameters.insertion_info.insertion_plane_normal_vector);
          find_centers_of_inplane_cells();
          mark_for_update = false;
        }

      MPI_Comm           communicator = triangulation.get_communicator();
      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);

      unsigned int this_mpi_process =
        Utilities::MPI::this_mpi_process(communicator);
      unsigned int total_mpi_process =
        Utilities::MPI::n_mpi_processes(communicator);
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

      // Processor zero knows how many particles every processor wants to insert
      auto number_of_particles_to_insert_per_core =
        Utilities::MPI::gather(communicator, empty_cells_on_proc.size(), 0);

      std::vector<int> starting_IDs_on_every_proc(total_mpi_process);
      std::vector<int> remained_particles_for_every_proc(total_mpi_process);
      // Only processor zero does this
      if (this_mpi_process == 0)
        {
          int starting_id_on_proc =
            particle_handler.get_next_free_particle_index();

          // we loop over the vector that we have created using the gather.
          // Every integer in that vector gives the number of empty_cell of a
          // certain processor.
          for (unsigned int i = 0; i < total_mpi_process; ++i)
            {
              // We subtract those integer from the following variable
              remained_particles_of_each_type -=
                number_of_particles_to_insert_per_core[i];

              // If the remained_particles_of_each_type variable become
              // negative, this means that we are trying to insert to many
              // particle in the simulation. Thus, we need to decrease the
              // number of empty_cell we want to insert in at this time step.
              if (remained_particles_of_each_type < 0)
                {
                  // At this point, we know the current processor have to many
                  // available cell to insert in (or just the right amount). The
                  // excess number of available cell on this processor is equal
                  // to the positive value of the
                  // remained_particles_of_each_type variable since that
                  // variable became negative on this processor.

                  // Thus, we decrease the number of particle to insert on the
                  // current processor. To do so, we add the negative variable
                  // to the number of available cells which guarantee that the
                  // number of inserted particles in the simulation will reach
                  // its maximum on this processor.

                  number_of_particles_to_insert_per_core[i] +=
                    remained_particles_of_each_type;

                  // We put the remained_particles_of_each_type variable  back
                  // to zero since we have just lower the number of cell
                  // available on this processor by the exceeding amount.
                  remained_particles_of_each_type = 0;

                  starting_IDs_on_every_proc[i] = starting_id_on_proc;
                  starting_id_on_proc +=
                    number_of_particles_to_insert_per_core[i];

                  // Loop over the remaining processor and put their number of
                  // available cell to zero.
                  for (unsigned int j = i + 1; j < (total_mpi_process); ++j)
                    {
                      number_of_particles_to_insert_per_core[j] = 0;
                      starting_IDs_on_every_proc[j]             = 0;
                    }
                  break;
                }
              starting_IDs_on_every_proc[i] = starting_id_on_proc;
              starting_id_on_proc += number_of_particles_to_insert_per_core[i];
            }
          std::fill(remained_particles_for_every_proc.begin(),
                    remained_particles_for_every_proc.end(),
                    remained_particles_of_each_type);
        }

      // Now, processor zero knows how many particles will be inserted by every
      // processor without exceeding the maximum number of particles.

      // We now need to scatter 3 information processor, first how
      // many particles it needs to insert. Second what is the ID of the
      // first particle it will insert at this insertion step. The starting ID
      // is mandatory, since we don't want two particles to have the same ID,
      // and we don't want to use the global_insertion method because it uses an
      // all_gather, which has a high computational cost. Third, the
      // remained_particles_of_each_type variable need to be updated proc other
      // than proc 0 will enter the insert function dans the simulation will
      // stop.

      // This scatters the number of particles to insert
      unsigned int number_of_particles_to_insert_on_this_core =
        Utilities::MPI::scatter(communicator,
                                number_of_particles_to_insert_per_core,
                                0);

      // This scatters the starting ID on every proc for the insertion
      unsigned int starting_ID_on_proc =
        Utilities::MPI::scatter(communicator, starting_IDs_on_every_proc, 0);

      // This scatters remained_particles_of_each_type on every proc
      remained_particles_of_each_type =
        Utilities::MPI::scatter(communicator,
                                remained_particles_for_every_proc,
                                0);

      // We yet didn't say which empty cells won't be use anymore.
      // empty_cells_on_proc needs to be shorter for some processor.

      // There's probably a better way without using a while loop, but we enter
      // in this while only once.
      while (empty_cells_on_proc.size() >
             number_of_particles_to_insert_on_this_core)
        {
          auto it =
            empty_cells_on_proc.end();   // Get an iterator to the last element
          std::advance(it, -1);          // Move the iterator back by one step
          empty_cells_on_proc.erase(it); // Erase the last element
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

      // Loop over the empty cells we have kept.
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
                                           starting_ID_on_proc++,
                                           cell,
                                           properties_of_one_particle);
        }
    }
}

template class PlaneInsertion<2>;

template class PlaneInsertion<3>;
