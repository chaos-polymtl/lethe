// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>
#include <core/lethe_grid_tools.h>

#include <dem/find_cell_neighbors.h>
#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/ray_tracing.h>
#include <dem/read_mesh.h>

#include <sys/stat.h>

template <int dim>
RayTracingSolver<dim>::RayTracingSolver(
  RayTracingSolverParameters<dim> &parameters,
  DEMSolverParameters<dim>        &dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , parameters(parameters)
  , dem_parameters(dem_parameters)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , particle_handler(triangulation,
                     mapping,
                     DEMProperties::PropertiesIndex::n_properties)
  , photon_handler(triangulation, mapping, 6)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , background_dh(triangulation)
{
  AssertThrow(!(parameters.model_parameters.load_balance_method ==
                  Parameters::Lagrangian::ModelParameters<
                    dim>::LoadBalanceMethod::dynamic_with_sparse_contacts ||
                parameters.model_parameters.load_balance_method ==
                  Parameters::Lagrangian::ModelParameters<
                    dim>::LoadBalanceMethod::dynamic),
              dealii::ExcMessage(
                "The \"dynamic\" and \"dynamic_with_sparse_contacts\""
                "load balancing methods are not supported"));
}

template <int dim>
void
RayTracingSolver<dim>::setup_parameters()
{
  // Print simulation starting information
  pcout << std::endl;
  std::stringstream ss;
  ss << "Running on " << n_mpi_processes << " rank(s)";
  announce_string(pcout, ss.str(), '*');

  // Check if the output directory exists
  std::string output_dir_name = parameters.simulation_control.output_folder;
  struct stat buffer;

  // If output directory does not exist, create it
  if (this_mpi_process == 0)
    {
      if (stat(output_dir_name.c_str(), &buffer) != 0)
        {
          create_output_folder(output_dir_name);
        }
    }

  // Get the pointer of the only instance of the action manager
  action_manager = DEMActionManager::get_action_manager();

  // Set the simulation control as SimulationControlAdjointSteady
  simulation_control =
    std::make_shared<SimulationControlSteady>(parameters.simulation_control);

  // Setup load balancing parameters and attach the correct functions to the
  // signals inside the triangulation
  // parameters.model_parameters.load_balance_method =
  // Parameters::Lagrangian::ModelParameters<dim>::frequent;
  load_balancing.set_parameters(parameters.model_parameters);

  AdaptiveSparseContacts<dim, DEMProperties::PropertiesIndex>
    dummy_sparse_contacts_object;
  load_balancing.copy_references(simulation_control,
                                 triangulation,
                                 photon_handler,
                                 dummy_sparse_contacts_object);
  load_balancing.connect_weight_signals();
}

template <int dim>
std::shared_ptr<Insertion<dim, DEMProperties::PropertiesIndex>>
RayTracingSolver<dim>::set_particle_insertion_type()
{
  using namespace Parameters::Lagrangian;
  typename InsertionInfo<dim>::InsertionMethod insertion_method =
    parameters.particle_insertion_info.insertion_method;

  std::vector<std::shared_ptr<Distribution>>
    dummy_size_distribution_object_container;

  switch (insertion_method)
    {
      case InsertionInfo<dim>::InsertionMethod::file:
        return std::make_shared<
          InsertionFile<dim, DEMProperties::PropertiesIndex>>(
          dummy_size_distribution_object_container,
          triangulation,
          dem_parameters);
      case InsertionInfo<dim>::InsertionMethod::list:
        return std::make_shared<
          InsertionList<dim, DEMProperties::PropertiesIndex>>(
          dummy_size_distribution_object_container,
          triangulation,
          dem_parameters);
      default:
        AssertThrow(
          false,
          dealii::ExcMessage(
            "For lethe-particles-ray-tracing simulation, only the list"
            " and file insertion are supported."));
        return nullptr;
    }
}


template <int dim>
void
RayTracingSolver<dim>::setup_background_dofs()
{
  FE_Q<dim> background_fe(1);
  background_dh.distribute_dofs(background_fe);
}

template <int dim>
void
RayTracingSolver<dim>::load_balance()
{
  load_balancing.check_load_balance_iteration();

  if (!action_manager->check_load_balance())
    return;

  TimerOutput::Scope t(this->computing_timer, "Load balancing");
  // Prepare particle handler for the adaptation of the triangulation to the
  // load
  photon_handler.prepare_for_coarsening_and_refinement();
  particle_handler.prepare_for_coarsening_and_refinement();

  pcout << "-->Repartitioning triangulation" << std::endl;
  triangulation.repartition();

  // Unpack the photon handler after the mesh has been repartitioned
  photon_handler.unpack_after_coarsening_and_refinement();
  particle_handler.unpack_after_coarsening_and_refinement();

  const auto average_minimum_maximum_cells =
    Utilities::MPI::min_max_avg(triangulation.n_active_cells(),
                                mpi_communicator);

  const auto average_minimum_maximum_photons =
    Utilities::MPI::min_max_avg(photon_handler.n_locally_owned_particles(),
                                mpi_communicator);

  pcout << "Load balance finished " << std::endl;
  pcout
    << "Average, minimum and maximum number of photons on the processors are "
    << average_minimum_maximum_photons.avg << " , "
    << average_minimum_maximum_photons.min << " and "
    << average_minimum_maximum_photons.max << std::endl;
  pcout << "Minimum and maximum number of cells owned by the processors are "
        << average_minimum_maximum_cells.min << " and "
        << average_minimum_maximum_cells.max << std::endl;

  setup_background_dofs();

  // Particle don't move, thus we sort them at the end of a load balance.
  // Photons are being sorted every pseudo time step, thus no need to sort them.
  particle_handler.sort_particles_into_subdomains_and_cells();

  // Exchange ghost particles
  particle_handler.exchange_ghost_particles(true);

  cells_local_neighbor_list.clear();
  cells_ghost_neighbor_list.clear();

  // Update cell neighbors
  find_cell_neighbors<dim, true>(triangulation,
                                 cells_local_neighbor_list,
                                 cells_ghost_neighbor_list);
}

template <int dim>
void
RayTracingSolver<dim>::print_insertion_info(const unsigned int &inserted_photon,
                                            const ConditionalOStream &pcout)
{
  std::stringstream ss;

  ss << inserted_photon << " photons "
     << " were inserted.";

  announce_string(pcout, ss.str(), '*');
}

template <int dim>
void
RayTracingSolver<dim>::insert_particles_and_photons()
{
  // Insert particles using the insertion object.
  particle_insertion_object->insert(particle_handler,
                                    triangulation,
                                    dem_parameters);

  // A vector of vectors, which contains all the properties of every photon.
  std::vector<Point<3>> insertion_points_on_proc;

  // For each photon, we store the initial location of the photon as a particle
  // property. This initial location is used to identify the right intersection
  // point when multiple intersections are found. The logic used therein is that
  // the intersection point which is closest to the initial location is the
  // correct one. We also store the displacement direction unit vector in this
  // vector since each photon has its own slightly off set from the
  std::vector<std::vector<double>> photon_properties;

  // Create variables for readability
  const std::vector<unsigned int> &n_photons_each_directions =
    parameters.ray_tracing_info.n_photons_each_directions;

  // Starting point
  const Point<dim> starting_insertion_point =
    parameters.ray_tracing_info.starting_point;

  // Insertion directions
  const Tensor<1, dim> first_dir =
    parameters.ray_tracing_info.insertion_directions_units_vector.at(0);
  const Tensor<1, dim> second_dir =
    parameters.ray_tracing_info.insertion_directions_units_vector.at(1);
  const Tensor<1, dim> third_dir =
    parameters.ray_tracing_info.insertion_directions_units_vector.at(2);

  // Steps
  const double step_first_dir =
    parameters.ray_tracing_info.step_between_photons_each_directions.at(0);
  const double step_second_dir =
    parameters.ray_tracing_info.step_between_photons_each_directions.at(1);
  const double step_third_dir =
    parameters.ray_tracing_info.step_between_photons_each_directions.at(2);

  // Reference displacement direction
  const Tensor<1, 3> ref_displacement_dir =
    parameters.ray_tracing_info.ref_displacement_tensor_unit;

  // Total number of photons to insert
  const unsigned int n_total_photons_to_insert =
    n_photons_each_directions.at(0) * n_photons_each_directions.at(1) *
    n_photons_each_directions.at(2);

  // Distributing photons between processors
  const unsigned int base_photons_per_proc =
    n_total_photons_to_insert / n_mpi_processes;
  const unsigned int remainder = n_total_photons_to_insert % n_mpi_processes;

  const unsigned int n_photons_to_insert_this_proc =
    base_photons_per_proc + (this_mpi_process < remainder ? 1 : 0);

  // First and last photon id on this processor. This will be used to find the
  // position of the photons.
  const unsigned int first_id =
    this_mpi_process * base_photons_per_proc +
    (this_mpi_process < remainder ? this_mpi_process : remainder);
  const unsigned int last_id = first_id + n_photons_to_insert_this_proc;

  // Generate the random offset in the insertion position and the
  // displacement direction
  std::vector<double> random_number_angular_1; // From 0 to 2*pi
  std::vector<double> random_number_angular_2; // From 0 to max_angular_offset
  std::vector<double> random_number_position;  // From 0 to max_insertion_offset

  random_number_angular_1.reserve(n_photons_to_insert_this_proc);
  random_number_angular_2.reserve(n_photons_to_insert_this_proc);
  random_number_position.reserve(3 * n_photons_to_insert_this_proc);

  create_random_number_container(
    random_number_angular_1,
    n_photons_to_insert_this_proc,
    2.0 * M_PI,
    parameters.ray_tracing_info.prn_seed_photon_displacement);

  create_random_number_container(
    random_number_angular_2,
    n_photons_to_insert_this_proc,
    parameters.ray_tracing_info.max_angular_offset,
    parameters.ray_tracing_info.prn_seed_photon_displacement);

  create_random_number_container(
    random_number_position,
    3 * n_photons_to_insert_this_proc,
    parameters.ray_tracing_info.max_insertion_offset,
    parameters.ray_tracing_info.prn_seed_photon_insertion);

  // Create the insertion location of every photon from their IDs.
  // Temporary variable
  Point<3>            temp_point;
  Tensor<1, 3>        temp_dir;
  std::vector<double> properties_of_one_photon(6);

  // Prepare the containers.
  insertion_points_on_proc.reserve(n_photons_to_insert_this_proc);
  photon_properties.resize(n_photons_to_insert_this_proc,
                           std::vector<double>(6));
  for (unsigned int id = first_id; id < last_id; ++id)
    {
      const unsigned int iz = id / (n_photons_each_directions.at(0) *
                                    n_photons_each_directions.at(1));

      const unsigned int rem = id % (n_photons_each_directions.at(0) *
                                     n_photons_each_directions.at(1));
      const unsigned int iy  = rem / n_photons_each_directions.at(0);
      const unsigned int ix  = rem % n_photons_each_directions.at(0);

      // ID relative to this processor
      const unsigned int id_on_proc = id - first_id;
      temp_point =
        starting_insertion_point +
        (ix * step_first_dir + random_number_position.at(3 * id_on_proc)) *
          first_dir +
        (iy * step_second_dir + random_number_position.at(3 * id_on_proc + 1)) *
          second_dir +
        (iz * step_third_dir + random_number_position.at(3 * id_on_proc + 2)) *
          third_dir;
      insertion_points_on_proc.push_back(temp_point);

      // For the displacement direction randomness, we need to find a vector
      // normal to ref_displacement_dir.
      const Tensor<1, 3> temp_normal_vector =
        (std::fabs(ref_displacement_dir[0]) < 0.9) ? Tensor<1, 3>({1, 0, 0}) :
                                                     Tensor<1, 3>({0, 1, 0});

      Tensor<1, 3> u =
        temp_normal_vector -
        temp_normal_vector *
          scalar_product(temp_normal_vector, ref_displacement_dir);
      u = u / u.norm();

      const Tensor<1, 3> v = cross_product_3d(ref_displacement_dir, u);
      temp_dir             = ref_displacement_dir +
                 random_number_angular_2.at(id_on_proc) *
                   (std::cos(random_number_angular_1.at(id_on_proc)) * u +
                    std::sin(random_number_angular_1.at(id_on_proc)) * v);

      for (unsigned int d = 0; d < 3; ++d)
        {
          properties_of_one_photon[d]     = temp_point[d];
          properties_of_one_photon[d + 3] = temp_dir[d];
        }
      photon_properties[id_on_proc] = properties_of_one_photon;
    }
  MPI_Comm communicator = triangulation.get_mpi_communicator();
  // Obtain global bounding boxes
  const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
    triangulation, IteratorFilters::LocallyOwnedCell());
  const auto global_bounding_boxes =
    Utilities::MPI::all_gather(communicator, my_bounding_box);

  // Insert the photons using the points and the assigned properties.
  photon_handler.insert_global_particles(insertion_points_on_proc,
                                         global_bounding_boxes,
                                         photon_properties);

  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);
  this->print_insertion_info(n_total_photons_to_insert, pcout);
}

template <int dim>
void
RayTracingSolver<dim>::write_output_results(
  const std::vector<Point<dim>> &points,
  const std::string             &folder,
  const std::string             &file_name)
{
  // Flatten local points into a buffer of chars (text format)
  std::ostringstream oss;
  if (this_mpi_process == 0)
    oss << "x, y, z\n";

  for (const auto &p : points)
    {
      oss << p[0] << "," << p[1] << "," << p[2] << "\n";
    }

  std::string local_str  = oss.str();
  int         local_size = local_str.size();

  // Gather sizes to compute offsets
  auto recv_counts = Utilities::MPI::all_gather(mpi_communicator, local_size);

  std::vector<int> off_set(n_mpi_processes, 0);
  int              total_size = 0;
  for (unsigned int i = 0; i < n_mpi_processes; ++i)
    {
      off_set[i] = total_size;
      total_size += recv_counts[i];
    }

  // Build full filename (same for all ranks)
  std::string full_filename = folder + "/" + file_name + ".xyz";

  // Open MPI file
  MPI_File fh;
  MPI_File_open(mpi_communicator,
                full_filename.c_str(),
                MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL,
                &fh);

  // Each rank writes at its own offset
  MPI_File_write_at(fh,
                    off_set[this_mpi_process],
                    local_str.data(),
                    (int)local_size,
                    MPI_CHAR,
                    MPI_STATUS_IGNORE);

  MPI_File_close(&fh);
}
template <int dim>
void
RayTracingSolver<dim>::finish_simulation(
  std::vector<Point<3>> &intersection_points)
{
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    this->computing_timer.print_summary();

  // Testing
  if (parameters.test.enabled)
    {
      // Proc 0 gathers all intersection points
      std::vector<std::vector<Point<3>>> all_intersection_points =
        Utilities::MPI::gather(mpi_communicator, intersection_points, 0);

      pcout << "x, y, z" << std::endl;
      // Loop over all intersection points
      for (const std::vector<Point<3>> &sub_vector : all_intersection_points)
        {
          for (const Point<3> &p : sub_vector)
            pcout << p << std::endl;
        }
    }
}

template <int dim>
template <bool move_photon>
void
RayTracingSolver<dim>::find_intersection(
  typename dem_data_structures<dim>::cells_neighbor_list &cell_list,
  ankerl::unordered_dense::map<
    types::particle_index,
    std::tuple<double, Point<dim>, Particles::ParticleIterator<dim>>>
    &photon_intersection_points_map)
{
  for (auto cell_neighbor_list_iterator = cell_list.begin();
       cell_neighbor_list_iterator != cell_list.end();
       ++cell_neighbor_list_iterator)
    {
      // The main cell
      const auto main_cell = cell_neighbor_list_iterator->begin();

      // Photons in main cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        photons_in_main_cell = photon_handler.particles_in_cell(*main_cell);

      // Loop on each photon in the main cell
      for (auto current_photon = photons_in_main_cell.begin();
           current_photon != photons_in_main_cell.end();
           ++current_photon)
        {
          // Get the photon properties. We need the properties to find the
          // insertion position of the current photon. This will be used to
          // differentiate between intersection points when there's more
          // than one.
          auto photon_properties = current_photon->get_properties();

          // Current photon insertion point.
          const Point<3> photon_insertion_point(
            {photon_properties[0], photon_properties[1], photon_properties[2]});

          // Displacement vector of the current photon
          Tensor<1, 3> photon_displacement_vector(
            {photon_properties[3], photon_properties[4], photon_properties[5]});

          const Point<dim> current_photon_location =
            current_photon->get_location();

          // Loop over the neighboring cells of the main cell. This includes
          // the main cell itself.
          auto &cell_neighbor_list = *cell_neighbor_list_iterator;
          for (auto current_neighboring_cell = cell_neighbor_list.begin();
               current_neighboring_cell != cell_neighbor_list.end();
               ++current_neighboring_cell)
            {
              // Particle iterators in the current neighboring cell
              typename Particles::ParticleHandler<dim>::particle_iterator_range
                particle_in_neighboring_cell =
                  particle_handler.particles_in_cell(*current_neighboring_cell);

              // This vector will be reused for every line-sphere
              // intersection search. There is maximum 2 intersection points
              // between a line and a sphere.
              std::vector<Point<dim>> current_intersection_points;
              current_intersection_points.reserve(2);

              // Loop of the particles in the current neighboring cells
              for (auto current_particle = particle_in_neighboring_cell.begin();
                   current_particle != particle_in_neighboring_cell.end();
                   ++current_particle)
                {
                  // Check if the current photon intersects the current
                  // particle
                  const Point<dim> particle_position =
                    current_particle->get_location();
                  const double particle_diameter =
                    current_particle
                      ->get_properties()[DEMProperties::PropertiesIndex::dp];

                  current_intersection_points =
                    LetheGridTools::find_line_sphere_intersection(
                      photon_insertion_point,
                      photon_displacement_vector,
                      particle_position,
                      particle_diameter);

                  // If the size of the vector is 0, there's no intersection
                  // point, thus we go to the next particle.
                  if (current_intersection_points.size() == 0)
                    continue;

                  // Otherwise, there is at least one intersection point.
                  // Initializing the needed variables.
                  double     new_distance;
                  Point<dim> new_closest_point;

                  // If the size of the vector is 1, we have one
                  // intersection point.
                  if (current_intersection_points.size() == 1)
                    {
                      new_distance = (photon_insertion_point -
                                      current_intersection_points[0])
                                       .norm();
                      new_closest_point = current_intersection_points[0];
                    }
                  // If the size is equal to 2, we need to find which one
                  // between the two new intersection point is the closest
                  // to the photon intersection point.
                  else
                    {
                      const double distance_0 = (photon_insertion_point -
                                                 current_intersection_points[0])
                                                  .norm();
                      const double distance_1 = (photon_insertion_point -
                                                 current_intersection_points[1])
                                                  .norm();

                      if (distance_0 < distance_1)
                        {
                          new_distance      = distance_0;
                          new_closest_point = current_intersection_points[0];
                        }
                      else
                        {
                          new_distance      = distance_1;
                          new_closest_point = current_intersection_points[1];
                        }
                    }
                  // We need to add the intersection point in the
                  // photon_insertion_point map. We first need to check
                  // if an intersection point already exist in the map
                  // for that given photon.
                  if (photon_intersection_points_map.find(
                        current_photon->get_id()) ==
                      photon_intersection_points_map.end())
                    {
                      // If it doesn't exist, we just add the new one.
                      photon_intersection_points_map[current_photon->get_id()] =
                        {new_distance, new_closest_point, current_photon};
                    }
                  else
                    {
                      // If it does exist, we need to check which
                      // intersection point is the closest to the
                      // current photon insertion point. We replace the
                      // value only if the new one is closer.
                      const double old_distance = std::get<0>(
                        photon_intersection_points_map[current_photon
                                                         ->get_id()]);

                      if (new_distance < old_distance)
                        {
                          // if the new distance is smaller, we replace the
                          // old one in the map. Otherwise we do nothing.
                          photon_intersection_points_map[current_photon
                                                           ->get_id()] = {
                            new_distance, new_closest_point, current_photon};
                        }
                    }
                }
            }
          // Even if we are removing photon at the end on this pseudo-time
          // step, we move the particle here since we are looping on each of
          // them in this loop. This way, we don't need to loop on each of
          // them at the end.
          if constexpr (move_photon)
            {
              const Point<dim> new_photon_location =
                current_photon_location +
                photon_displacement_vector * displacement_distance;
              current_photon->set_location(new_photon_location);
              // std::cout<< new_photon_location << std::endl;
            }
        }
    }
}


template <int dim>
void
RayTracingSolver<dim>::solve()
{
  // Set up the parameters
  setup_parameters();

  // Reading mesh
  read_mesh(parameters.mesh,
            action_manager->check_restart_simulation(),
            pcout,
            triangulation,
            dem_parameters.boundary_conditions);

  displacement_distance = 0.5 * GridTools::minimal_cell_diameter(triangulation);

  // Set particle insertion object
  particle_insertion_object = set_particle_insertion_type();

  // Insert photon and particles
  insert_particles_and_photons();

  // Particles won't move, thus we can sort them.
  particle_handler.sort_particles_into_subdomains_and_cells();
  // Exchange ghost particles
  particle_handler.exchange_ghost_particles(true);

  // Every intersection point from the beginning of the simulation.
  std::vector<Point<dim>> total_intersection_points;
  total_intersection_points.reserve(photon_handler.n_global_particles());

  // Map :
  // <photon_iterator, (distance, intersection point<dim>, iterator to remove)>
  ankerl::unordered_dense::map<
    types::particle_index,
    std::tuple<double, Point<dim>, Particles::ParticleIterator<dim>>>
    photon_intersection_points_map;

  photon_intersection_points_map.reserve(
    photon_handler.n_locally_owned_particles());

  // Update cell neighbors
  find_cell_neighbors<dim, true>(triangulation,
                                 cells_local_neighbor_list,
                                 cells_ghost_neighbor_list);

  // Particle don't move, thus we sort them at the end of a load balance.
  // Photons are being sorted every pseudo time step, thus no need to sort them.
  particle_handler.sort_particles_into_subdomains_and_cells();

  // Exchange ghost particles
  particle_handler.exchange_ghost_particles(true);

  while (photon_handler.n_global_particles() != 0)
    {
      simulation_control->increment_iteration();

      pcout << "Remaining photon : " << photon_handler.n_global_particles()
            << std::endl;

      // Load balancing (if needed)
      load_balance();
      // Since last pseudo time step, there is a high chance that the photons
      // have changed cell considering the size of their displacement relative
      // to the cell size. Thus, we sort them.

      photon_handler.sort_particles_into_subdomains_and_cells();

      // Loop over each local cell in the triangulation. To do this, we loop
      // over each cell neighbor list. Each local cell has a list and the first
      // iterator in the main cell itself.
      find_intersection<true>(cells_local_neighbor_list,
                              photon_intersection_points_map);

      find_intersection<false>(cells_ghost_neighbor_list,
                               photon_intersection_points_map);

      // Remove all the photon that have found their intersection. In other
      // words, every photon that is present in the map. We also store the
      // intersection points in the appropriate vector.
      std::vector<Particles::ParticleIterator<dim>> photon_iterators_to_remove;
      photon_iterators_to_remove.reserve(photon_intersection_points_map.size());

      for (auto it = photon_intersection_points_map.begin();
           it != photon_intersection_points_map.end();
           ++it)
        {
          auto photon_tuple = it->second;
          total_intersection_points.push_back(std::get<1>(photon_tuple));
          photon_iterators_to_remove.push_back(std::get<2>(photon_tuple));
        }

      // Removes the photon and clear the containers.
      photon_handler.remove_particles(photon_iterators_to_remove);
      photon_iterators_to_remove.clear();
      photon_intersection_points_map.clear();

      action_manager->reset_triggers();
    }
  // std::cout << this_mpi_process << " " << total_intersection_points.size()
  //           << std::endl;
  write_output_results(total_intersection_points,
                       parameters.simulation_control.output_folder,
                       parameters.simulation_control.output_name);

  finish_simulation(total_intersection_points);
}


// template class RayTracingSolver<2>;
template class RayTracingSolver<3>;
