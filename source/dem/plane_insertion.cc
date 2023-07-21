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
  const parallel::distributed::Triangulation<dim> &triangulation)
  : remained_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
{
  // Initializing current inserting particle type
  current_inserting_particle_type = 0;
  this->inserted_this_step        = 0;
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

          // Vector that goes from the plane reference point to the vertice 0.
          Tensor<1, 3> connecting_vector_ref =
            point_nd_to_3d(cell->vertex(0)) - plane_point;

          // Normal distance between vertex 0 and the plan
          double vertex_wall_distance_ref =
            std::abs(connecting_vector_ref * plane_normal_vector);

          // Loop over all the vertices of the cell
          for (unsigned int vertex_id = 1; vertex_id < cell->n_vertices();
               ++vertex_id)
            {
              // Vector that goes from the plane reference point to the vertice
              // n-th
              Tensor<1, 3> connecting_vector =
                point_nd_to_3d(cell->vertex(vertex_id)) - plane_point;

              // Normal distance between vertex n-th and the plan
              double vertex_wall_distance =
                std::abs(connecting_vector * plane_normal_vector);

              // If the multiplication of the two distances is negatif, then at
              // least one of vertex is on both sides of the plane. If the
              // multiplication is equal to zero, then one vertex is on the
              // plane.
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
      current_inserting_particle_type !=
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
      auto n_mpi_process    = Utilities::MPI::n_mpi_processes(communicator);

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // List of cells to insert at this time step
      std::vector<Point<dim>> insertion_points_on_proc_this_step;

      // Loop over the cell inplane
      for (const auto &cell : plane_cells_for_insertion)
        {
          // if the cell is empty and is locally own
          if (particle_handler.n_particles_in_cell(cell) == 0)
            {
              // Position of the center of the cell (To remove
              Point<dim> particle_location =
                cells_centers.at(cell->active_cell_index());

              // Generate the random point of insertion
              particle_location(0) +=
                static_cast<double>(rand()) / static_cast<double>(RAND_MAX) *
                dem_parameters.insertion_info.random_number_range;

              particle_location(1) +=
                static_cast<double>(rand()) / static_cast<double>(RAND_MAX) *
                dem_parameters.insertion_info.random_number_range;

              if constexpr (dim == 3)
                {
                  particle_location(2) +=
                    static_cast<double>(rand()) /
                    static_cast<double>(RAND_MAX) *
                    dem_parameters.insertion_info.random_number_range;
                }
              // Point to insert at this time step
              insertion_points_on_proc_this_step.push_back(particle_location);
            }
        }
      // Gather to core 0 the number of particle every core can insert

      auto number_of_particles_to_insert_per_core =
        Utilities::MPI::all_gather(communicator,
                                   insertion_points_on_proc_this_step.size());

      for (unsigned int i = 0;
           i < number_of_particles_to_insert_per_core.size();
           ++i)
        {
          remained_particles_of_each_type -=
            number_of_particles_to_insert_per_core[i];
          // If the remaining number of particle become negative, this means
          // that we want to insert to many particle at this time step. So...

          if (remained_particles_of_each_type <= 0)
            {
              // ... we decrease the number of particle at the current core in
              // the sum and...


              number_of_particles_to_insert_per_core[i] +=
                remained_particles_of_each_type;


              remained_particles_of_each_type = 0;

              for (unsigned int j = i + 1;
                   j < number_of_particles_to_insert_per_core.size();
                   ++j)
                {
                  // we put the remaining cores to zero
                  number_of_particles_to_insert_per_core[j] = 0;
                }
              break;
            }
        }

      std::vector<Point<dim>> new_insertion_points_on_proc(
        insertion_points_on_proc_this_step.begin(),
        insertion_points_on_proc_this_step.begin() +
          number_of_particles_to_insert_per_core[this_mpi_process]);



      this->assign_particle_properties(dem_parameters,
                                       new_insertion_points_on_proc.size(),
                                       current_inserting_particle_type,
                                       this->particle_properties);

      particle_handler.insert_global_particles(new_insertion_points_on_proc,
                                               global_bounding_boxes,
                                               this->particle_properties);

      //      particle_handler.insert_particle(particle_location,
      //                                       ref_point,
      //                                       particle_counter++,
      //                                       cell,
      //                                       properties_of_one_particle);
    }
}

;
template class PlaneInsertion<2>;
template class PlaneInsertion<3>;
