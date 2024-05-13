#include <core/utilities.h>

#include <dem/insertion_clear_and_add.h>

using namespace DEM;

template <int dim>
InsertionClearAndAdd<dim>::InsertionClearAndAdd(
  const DEMSolverParameters<dim>                  &dem_parameters,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const std::vector<std::shared_ptr<Distribution>>
    &distribution_object_container)
  : Insertion<dim>(distribution_object_container)
  , remaining_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
  , number_of_files(dem_parameters.insertion_info.list_of_input_files.size())
  , insertion_files(dem_parameters.insertion_info.list_of_input_files)
{
  // Boost signal for load balancing
  mark_for_update = true;
  this->change_to_triangulation =
    triangulation.signals.any_change.connect([&] { mark_for_update = true; });

  // Initializing current inserting particle type and file id
  this->current_inserting_particle_type = 0;
  this->current_file_id                 = 0;

  for (unsigned int i = 0; i < 3; ++i)
    {
      if (dem_parameters.insertion_info.clear_box_point_1[i] <=
          dem_parameters.insertion_info.clear_box_point_2[i])
        {
          p_min[i] = dem_parameters.insertion_info.clear_box_point_1[i];
          p_max[i] = dem_parameters.insertion_info.clear_box_point_2[i];
        }
      else
        {
          p_min[i] = dem_parameters.insertion_info.clear_box_point_2[i];
          p_max[i] = dem_parameters.insertion_info.clear_box_point_1[i];
        }
    }
}

template <int dim>
void
InsertionClearAndAdd<dim>::find_in_clearing_box_cells(
  const parallel::distributed::Triangulation<dim> &triangulation)
{
  // Clearing the containers
  in_the_clearing_box.clear();
  edge_of_box.clear();
  bool partially_inside, completely_inside;
  // Looping through cells
  for (const auto &cell : triangulation.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          partially_inside  = false;
          completely_inside = true;
          // Loop through its vertices
          for (unsigned int vertex_id = 0; vertex_id < cell->n_vertices();
               ++vertex_id)
            {
              Point<3> vertex = point_nd_to_3d(cell->vertex(vertex_id));

              // Check if the n_th vertex is in the box
              if (p_min[0] <= vertex[0] && vertex[0] <= p_max[0] &&
                  p_min[1] <= vertex[1] && vertex[1] <= p_max[1] &&
                  p_min[2] <= vertex[2] && vertex[2] <= p_max[2])
                {
                  partially_inside = true;
                }
              else
                {
                  completely_inside = false;
                }
            }
          // If the cell is completely inside the clearing box, both bool will
          // be true. We need to add the cell iterator to the
          // in_the_clearing_box container.

          // If the cell is partially inside, at least one of the "if" will have
          // failed, thus the second bool will be false but the first one true.
          // We need to add the cell iterator to the edge_of_box container.

          // If the cell is outside, both bool will be false, and we do nothing.
          if (completely_inside)
            {
              in_the_clearing_box.insert(cell);
            }
          else if (!completely_inside && partially_inside)
            {
              edge_of_box.insert(cell);
            }
        }
    }
}

template <int dim>
void
InsertionClearAndAdd<dim>::insert(
  Particles::ParticleHandler<dim>                 &particle_handler,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
{
  if (remaining_particles_of_each_type == 0 &&
      this->current_inserting_particle_type !=
        dem_parameters.lagrangian_physical_properties.particle_type_number - 1)
    {
      remaining_particles_of_each_type =
        dem_parameters.lagrangian_physical_properties.number.at(
          ++this->current_inserting_particle_type);
    }


  if (remaining_particles_of_each_type > 0)
    {
      // Check if the triangulation has changed
      if (mark_for_update)
        {
          find_in_clearing_box_cells(triangulation);
          mark_for_update = false;
        }

      // Vector containing every particle iterator to remove
      std::vector<
        typename dealii::Particles::ParticleHandler<dim>::particle_iterator>
        to_remove_iterators;

      // Reserve to the maximum number of particle on this proc
      to_remove_iterators.reserve(particle_handler.n_locally_owned_particles());

      // Loop over the first container
      for (const auto &cell_in_box : in_the_clearing_box)
        {
          // Check if this cell has particles
          auto particles_in_cell =
            particle_handler.particles_in_cell(cell_in_box);

          if (!particles_in_cell.empty())
            {
              // Loop over the particle in the cell
              for (auto particle_in_cell = particles_in_cell.begin();
                   particle_in_cell != particles_in_cell.end();
                   ++particle_in_cell)
                {
                  // Since we know the cell is fully inside the box, we can
                  // delete every particles in it.

                  to_remove_iterators.push_back(particle_in_cell);
                }
            }
        }

      // Loop over the second container
      for (auto cell_edge_of_box = edge_of_box.begin();
           cell_edge_of_box != edge_of_box.end();
           ++cell_edge_of_box)
        {
          // Check if this cell has particle
          auto particles_in_cell =
            particle_handler.particles_in_cell(*cell_edge_of_box);
          const bool particles_exist_in_cell = !particles_in_cell.empty();

          if (particles_exist_in_cell)
            {
              // Loop over the particles in the cell
              for (auto particle_in_cell = particles_in_cell.begin();
                   particle_in_cell != particles_in_cell.end();
                   ++particle_in_cell)
                {
                  // We need to check if the particle is inside since the cell
                  // is at the edge.
                  Point<3> particle_position =
                    point_nd_to_3d(particle_in_cell->get_location());
                  if (p_min[0] <= particle_position[0] &&
                      p_max[0] >= particle_position[0] &&
                      p_min[1] <= particle_position[1] &&
                      p_max[1] >= particle_position[1] &&
                      p_min[2] <= particle_position[2] &&
                      p_max[2] >= particle_position[2])
                    {
                      to_remove_iterators.push_back(particle_in_cell);
                    }
                }
            }
        }
      particle_handler.remove_particles(to_remove_iterators);

      // Read the input file
      std::map<std::string, std::vector<double>> particles_data;
      fill_vectors_from_file(particles_data,
                             insertion_files.at(current_file_id),
                             ";");

      // Number of particles in the file
      unsigned int n_total_particles_to_insert = particles_data["p_x"].size();

      // Adjusting the value in case we exceed the maximum number of particle in
      // the simulation.
      n_total_particles_to_insert =
        std::min(remaining_particles_of_each_type, n_total_particles_to_insert);

      // Processor 0 will be the only one inserting particles
      MPI_Comm communicator = triangulation.get_communicator();
      auto this_mpi_process = Utilities::MPI::this_mpi_process(communicator);
      const unsigned int n_particles_to_insert_this_proc =
        this_mpi_process == 0 ? n_total_particles_to_insert : 0;

      std::vector<Point<dim>> insertion_points_on_proc_this_step;
      insertion_points_on_proc_this_step.reserve(
        n_particles_to_insert_this_proc);

      if (this_mpi_process == 0)
        {
          for (unsigned int p = 0; p < n_particles_to_insert_this_proc; ++p)
            {
              if constexpr (dim == 2)
                {
                  insertion_points_on_proc_this_step.emplace_back(Point<dim>(
                    {particles_data["p_x"][p], particles_data["p_y"][p]}));
                }

              if constexpr (dim == 3)
                {
                  insertion_points_on_proc_this_step.emplace_back(
                    Point<dim>({particles_data["p_x"][p],
                                particles_data["p_y"][p],
                                particles_data["p_z"][p]}));
                }
            }
        }

      // Obtain global bounding boxes
      const auto my_bounding_box =
        GridTools::compute_mesh_predicate_bounding_box(
          triangulation, IteratorFilters::LocallyOwnedCell());
      const auto global_bounding_boxes =
        Utilities::MPI::all_gather(communicator, my_bounding_box);

      // A vector of vectors, which contains all the properties of all particles
      // about to get inserted
      std::vector<std::vector<double>> particle_properties;

      // Assign inserted particles properties
      this->assign_particle_properties_for_file_insertion(
        dem_parameters,
        n_particles_to_insert_this_proc,
        particles_data,
        particle_properties);

      // Insert the particles using the points and assigned properties
      particle_handler.insert_global_particles(
        insertion_points_on_proc_this_step,
        global_bounding_boxes,
        particle_properties);

      // Update number of particle remaining to be inserted
      remaining_particles_of_each_type -= n_total_particles_to_insert;

      ConditionalOStream pcout(
        std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0);
      this->print_insertion_info(n_total_particles_to_insert,
                                 remaining_particles_of_each_type,
                                 this->current_inserting_particle_type,
                                 pcout);

      // Adjusting the file id
      current_file_id += 1;
      current_file_id = current_file_id % number_of_files;
    }
}

template <int dim>
void
InsertionClearAndAdd<dim>::assign_particle_properties_for_file_insertion(
  const DEMSolverParameters<dim>             &dem_parameters,
  const unsigned int                         &inserted_this_step_this_proc,
  std::map<std::string, std::vector<double>> &particles_data,
  std::vector<std::vector<double>>           &particle_properties)
{
  // Clearing and resizing particle_properties
  particle_properties.reserve(inserted_this_step_this_proc);

  // Getting properties as local parameters
  auto physical_properties = dem_parameters.lagrangian_physical_properties;

  // A loop is defined over the number of particles which are going to be
  // inserted at this step
  for (unsigned int particle_counter = 0;
       particle_counter < inserted_this_step_this_proc;
       ++particle_counter)
    {
      double type     = this->current_inserting_particle_type;
      double diameter = particles_data["diameters"][particle_counter];
      double density =
        physical_properties
          .density_particle[this->current_inserting_particle_type];
      double vel_x        = particles_data["v_x"][particle_counter];
      double vel_y        = particles_data["v_y"][particle_counter];
      double vel_z        = particles_data["v_z"][particle_counter];
      double omega_x      = particles_data["w_x"][particle_counter];
      double omega_y      = particles_data["w_y"][particle_counter];
      double omega_z      = particles_data["w_z"][particle_counter];
      double fem_force_x  = particles_data["fem_force_x"][particle_counter];
      double fem_force_y  = particles_data["fem_force_y"][particle_counter];
      double fem_force_z  = particles_data["fem_force_z"][particle_counter];
      double fem_torque_x = particles_data["fem_torque_x"][particle_counter];
      double fem_torque_y = particles_data["fem_torque_y"][particle_counter];
      double fem_torque_z = particles_data["fem_torque_z"][particle_counter];
      double mass         = density * 4. / 3. * M_PI *
                    Utilities::fixed_power<3, double>(diameter * 0.5);
      double volumetric_contribution = 0;

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

      particle_properties.push_back(properties_of_one_particle);
      properties_of_one_particle.clear();
    }
}
template class InsertionClearAndAdd<2>;
template class InsertionClearAndAdd<3>;
