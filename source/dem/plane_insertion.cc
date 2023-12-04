#include <core/tensors_and_points_dimension_manipulation.h>

#include <dem/plane_insertion.h>

using namespace DEM;

// The constructor of plane insertion class. In the constructor, we find which
// cells are going to be use for the insertion and we also find the centers of
// those cells.
template <int dim>
PlaneInsertion<dim>::PlaneInsertion(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<std::shared_ptr<Distribution>>
    &distribution_object_container)
  : Insertion<dim>(distribution_object_container)
  , particles_of_each_type_remaining(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;

  // Finding which cells are inplane
  this->find_inplane_cells(
    triangulation,
    dem_parameters.insertion_info.insertion_plane_point,
    dem_parameters.insertion_info.insertion_plane_normal_vector);

  // Finding the center of those cells
  this->find_centers_of_inplane_cells();
  mark_for_update = false;

  // Boost signal for load balancing
  change_to_triangulation =
    triangulation.signals.any_change.connect([&] { mark_for_update = true; });

  // Initializing the variable for the random position of particle
  maximum_range_for_randomness =
    dem_parameters.insertion_info.random_number_range /
    static_cast<double>(RAND_MAX);
}

template <int dim>
void
PlaneInsertion<dim>::find_inplane_cells(
  const parallel::distributed::Triangulation<dim> &triangulation,
  Point<3>                                         plane_point,
  Tensor<1, 3>                                     plane_normal_vector)
{
  plane_cells_for_insertion.clear();

  // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      // If the cell is owned by the processor
      if (cell->is_locally_owned())
        {
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

              // Normal distance between vertex n-th and the plane
              double vertex_wall_distance =
                connecting_vector * plane_normal_vector;

              // If the multiplication of the two distances is negative, then at
              // least one of vertex of the cell is on both sides of the plane.
              // If the multiplication is equal to zero, then one vertex is on
              // the plane.
              if (vertex_wall_distance_ref * vertex_wall_distance <= 0)
                {
                  plane_cells_for_insertion.insert(cell);
                  break; // If so, we break the loop, and go to the next cell
                }
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

// The main insertion function, insert_particle, is utilized to insert the
// particle at the cell center with a random shifts particles
template <int dim>
void
PlaneInsertion<dim>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  if (particles_of_each_type_remaining == 0 &&
      this->current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      particles_of_each_type_remaining =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++current_inserting_particle_type);
    }
  if (particles_of_each_type_remaining > 0)
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
              particles_of_each_type_remaining -=
                number_of_particles_to_insert_per_core[i];

              // If the particles_of_each_type_remaining variable become
              // negative, this means that we are trying to insert too many
              // particle in the simulation. Thus, we need to decrease the
              // number of empty_cell we want to insert in at this time step.
              if (particles_of_each_type_remaining < 0)
                {
                  // At this point, we know the current processor have to many
                  // available cell to insert in (or just the right amount). The
                  // excess number of available cell on this processor is equal
                  // to the positive value of the
                  // particles_of_each_type_remaining variable since that
                  // variable became negative on this processor.

                  // Thus, we decrease the number of particle to insert on the
                  // current processor. To do so, we add the negative variable
                  // to the number of available cells which guarantee that the
                  // number of inserted particles in the simulation will reach
                  // its maximum on this processor.

                  number_of_particles_to_insert_per_core[i] +=
                    particles_of_each_type_remaining;

                  // We put the particles_of_each_type_remaining variable  back
                  // to zero since we have just lower the number of cell
                  // available on this processor by the exceeding amount.
                  particles_of_each_type_remaining = 0;

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
                    particles_of_each_type_remaining);
        }

      // Now, processor zero knows how many particles will be inserted by every
      // processor without exceeding the maximum number of particles.

      // We now need to scatter 3 information per processor, first how
      // many particles it needs to insert. Second what is the ID of the
      // first particle it will insert at this insertion step. The starting ID
      // is mandatory, since we don't want two particles to have the same ID,
      // and we don't want to use the global_insertion method because it uses an
      // all_gather, which has a high computational cost. Third, the
      // particles_of_each_type_remaining variable needs to be updated on
      // processors other than processor 0 otherwise the insert function will
      // become incoherent.

      // This scatters the number of particles to insert
      unsigned int number_of_particles_to_insert_on_this_core =
        Utilities::MPI::scatter(communicator,
                                number_of_particles_to_insert_per_core,
                                0);

      // This scatters the starting ID on every proc for the insertion
      unsigned int starting_ID_on_proc =
        Utilities::MPI::scatter(communicator, starting_IDs_on_every_proc, 0);

      // This scatters particles_of_each_type_remaining on every proc
      particles_of_each_type_remaining =
        Utilities::MPI::scatter(communicator,
                                remained_particles_for_every_proc,
                                0);

      // We yet didn't choose which empty cells won't be used anymore.

      while (empty_cells_on_proc.size() >
             number_of_particles_to_insert_on_this_core)
        {
          auto it = empty_cells_on_proc
                      .begin(); // Get the iterator of the first element
          empty_cells_on_proc.erase(it); // Erase the first element
        }

      // A vector of vectors, which contains all the properties of all inserted
      // particles at each insertion step
      std::vector<std::vector<double>> particle_properties;

      this->assign_particle_properties(
        dem_parameters,
        number_of_particles_to_insert_on_this_core,
        current_inserting_particle_type,
        particle_properties);

      // This is to iterate over the particle_properties vector
      unsigned int i = 0;

      // Loop over the empty cells we have kept.
      for (const auto &cell : empty_cells_on_proc)
        {
          // Position of the center of the cell
          Point<dim> insertion_location =
            cells_centers.at(cell->active_cell_index());

          // Generate the random point of insertion
          insertion_location(0) +=
            static_cast<double>(rand()) * maximum_range_for_randomness;

          insertion_location(1) +=
            static_cast<double>(rand()) * maximum_range_for_randomness;

          if constexpr (dim == 3)
            {
              insertion_location(2) +=
                static_cast<double>(rand()) * maximum_range_for_randomness;
            }

          // Insertion
          Point<dim> ref_point;
          particle_handler.insert_particle(insertion_location,
                                           ref_point,
                                           starting_ID_on_proc++,
                                           cell,
                                           particle_properties.at(i++));
        }
    }
}

template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
