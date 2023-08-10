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
  particle_counter                = 0;

  // Finding which cells are inplane
  this->find_inplane_cells(
    triangulation,
    dem_parameters.insertion_info.insertion_plane_point,
    dem_parameters.insertion_info.insertion_plane_normal_vector);

  // // Finding the center of those cells
  this->find_centers_of_inplane_cells();
}

// This function creates a vector of random doubles using the input paramteres
// in the parameter handler
template <int dim>
void
PlaneInsertion<dim>::create_random_offset_container(
  std::vector<Point<dim>> &random_container,
  const double             random_number_range)
{
  for (unsigned int i = 0; i < this->inserted_this_step; ++i)
    {
      Point<dim> insertion_off_set;
      for (unsigned int j = 0; j < dim; ++j)
        {
          insertion_offset(j) += static_cast<double>(rand()) /
                                 static_cast<double>(RAND_MAX) *
                                 random_number_range;
        }
      random_container.emplace_back(insertion_off_set);
    }
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

      // Processor zero knows how many particles every processor wants to insert
      auto number_of_particles_to_insert_per_core =
        Utilities::MPI::gather(communicator, empty_cells_on_proc.size(), 0);

      std::vector<int> starting_IDs_on_every_proc;
      // Only processor zero does this
      if (this_mpi_process == 0)
        {
          int starting_id_on_proc =
            particle_handler.get_next_free_particle_index();

          // we loop over the vector that we have created using the gather.
          // Every integer in that vector gives the number of empty_cell of a
          // certain processor.
          for (unsigned int i = 0;
               i < number_of_particles_to_insert_per_core.size();
               ++i)
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

                  // Loop over the remaining processor and put their number of
                  // available cell to zero.
                  for (unsigned int j = i + 1;
                       j < number_of_particles_to_insert_per_core.size();
                       ++j)
                    {
                      number_of_particles_to_insert_per_core[j] = 0;
                    }
                  break;
                }
              starting_IDs_on_every_proc.emplace_back(starting_id_on_proc);
              starting_id_on_proc += number_of_particles_to_insert_per_core[i];
            }
        }
      // Now, processor zero knows how many particles will be inserted by every
      // processor without exceeding the maximum number of particles.

      // We now need to scatter two details to every processor, first how
      // many particles it needs to insert and second what is the ID of the
      // first particle it will insert at this insertion step.

      // The transmission of the starting ID is mandatory, since we don't
      // want two particles to have the same ID, and we don't want to use the
      // global_insertion method because it uses an all_gather, which has a
      // high computational cost.

      // This scatters the number of particles to insert
      int number_of_particles_to_insert_on_this_core =
        Utilities::MPI::scatter(communicator,
                                number_of_particles_to_insert_per_core,
                                0);

      // This scatters the starting ID on every proc for the insertion
      int starting_ID_on_proc =
        Utilities::MPI::scatter(communicator, starting_IDs_on_every_proc, 0);

      // Every processor knows how many particles it needs to insert at this
      // insertion this time step. Next step is to generate a vector Point<dim>
      // and a list of IDs to insert particle using the insert_particles
      // function.

      // Created the vectors of IDs with the right size
      std::vector<int> vector_IDs(number_of_particles_to_insert_on_this_core);
      std::vector<Point<dim>> vector_insertion_point(
        number_of_particles_to_insert_on_this_core);
      // For the IDs, goes from starting_ID_on_proc to starting_ID_on_proc +
      // number_of_particles_to_insert_on_this_core
      std::iota(vector_IDs.begin(), vector_IDs.end(), starting_ID_on_proc);

      // Create the offset Point<dim> vector
      std::vector<Point<dim>> random_insert_offset_vector;
      random_number_vector.reserve(number_of_particles_to_insert_on_this_core);

      // Call random Point generator
      this->create_random_offset_container(
        random_insert_offset_vector,
        dem_parameters.insertion_info.random_number_range);

      // Create the final insertion location vector
      std::vector<Point<dim>> insertion_location;
      insertion_location.reserve(number_of_particles_to_insert_on_this_core);

      while (empty_cells_on_proc.size() <
             number_of_particles_to_insert_per_core[this_mpi_process])
        {
          auto it = empty_cells_on_proc.end();
          empty_cells_on_proc.erase(it);
        }

      int i = 0;
      for (const auto &cell : empty_cells_on_proc)
        {
          insertion_location[i] = cells_centers.at(cell->active_cell_index()) +
                                  random_insert_offset_vector[i++];
        }

      particle_handler.insert_particles( insertion_location);
      // ICI OLIVIER!!!
    }
}
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
