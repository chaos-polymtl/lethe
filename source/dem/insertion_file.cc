#include <core/utilities.h>

#include <dem/insertion_file.h>


using namespace DEM;

template <int dim>
InsertionFile<dim>::InsertionFile(
  const std::vector<std::shared_ptr<Distribution>>
    &distribution_object_container,
  const parallel::distributed::Triangulation<dim> &triangulation,
  const DEMSolverParameters<dim>                  &dem_parameters)
  : Insertion<dim>(distribution_object_container, triangulation, dem_parameters)
  , remaining_particles_of_each_type(
      dem_parameters.lagrangian_physical_properties.number.at(0))
  , number_of_files(dem_parameters.insertion_info.list_of_input_files.size())
  , insertion_files(dem_parameters.insertion_info.list_of_input_files)
{
  // Initializing current inserting particle type
  this->current_inserting_particle_type = 0;

  // Initializing current inserting particle type and file id
  current_inserting_particle_type = 0;
  current_file_id                 = 0;
}
template <int dim>
void
InsertionFile<dim>::insert(
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

      if (this->removing_particles_in_region)
        {
          if (this->mark_for_update)
            {
              this->find_in_removing_box_cells(triangulation);
              this->mark_for_update = false;
            }
          this->remove_particles_in_box(particle_handler);
        }

      // Read the input file
      std::map<std::string, std::vector<double>> particles_data;
      fill_vectors_from_file(particles_data,
                             insertion_files.at(current_file_id),
                             ";");
      current_file_id++;
      current_file_id = current_file_id % number_of_files;

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
    }
}

template <int dim>
void
InsertionFile<dim>::assign_particle_properties_for_file_insertion(
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

template class InsertionFile<2>;
template class InsertionFile<3>;
