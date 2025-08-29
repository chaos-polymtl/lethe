// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

// core
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
  RayTracingSolverParameters<dim> parameters,
  DEMSolverParameters<dim>        dem_parameters)
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
  , photon_handler(triangulation, mapping, dim)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , background_dh(triangulation)
  , photon_displacement_vector(
      parameters.ray_tracing_info.photon_displacement_vector)
{
  if (parameters.model_parameters.load_balance_method ==
        Parameters::Lagrangian::ModelParameters<
          dim>::LoadBalanceMethod::dynamic_with_sparse_contacts ||
      parameters.model_parameters.load_balance_method ==
        Parameters::Lagrangian::ModelParameters<
          dim>::LoadBalanceMethod::dynamic)
    throw std::runtime_error(
      "The \"dynamic\" and \"dynamic_with_sparse_contacts\""
      "load balancing methods are not supported");
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
  simulation_control = std::make_shared<SimulationControlSteady>(
    parameters.simulation_control);

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
        {
          return std::make_shared<
            InsertionFile<dim, DEMProperties::PropertiesIndex>>(
            dummy_size_distribution_object_container,
            triangulation,
            dem_parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::list:
        {
          return std::make_shared<
            InsertionList<dim, DEMProperties::PropertiesIndex>>(
            dummy_size_distribution_object_container,
            triangulation,
            dem_parameters);
        }
      default:
        throw(std::runtime_error("Invalid insertion method."));
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

  pcout << "-->Repartitionning triangulation" << std::endl;
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

  // A vector which contains all the insertion point of every photon.
  std::vector<Point<dim>> insertion_points_on_proc;

  // A vector of vectors, which contains all the properties of all particles
  // about to get inserted.
  // Each photon has as its properties its insertion location. This way,
  // when many intersection points will be found during a pseudo timestep,
  // its properties will be used to identify the closest point from its
  // insertion insertion point.
  std::vector<std::vector<double>> photon_properties;

  unsigned int max_n_photon_first_dir =
    parameters.ray_tracing_info.number_of_photon_first_direction;

  unsigned int max_n_photon_second_dir =
    parameters.ray_tracing_info.number_of_photon_second_direction;

  // Processor 0 will be the only one inserting photon
  const unsigned int n_particles_to_insert_this_proc =
    this_mpi_process == 0 ? max_n_photon_first_dir * max_n_photon_second_dir :
                            0;

  insertion_points_on_proc.reserve(n_particles_to_insert_this_proc);

  photon_properties.reserve(n_particles_to_insert_this_proc);

  // Create the photon insertion location
  if (this_mpi_process == 0)
    {
      // Create variable for readability
      Point<dim> starting_insertion_point =
        parameters.ray_tracing_info.starting_point;
      Tensor<1, dim> first_dir  = parameters.ray_tracing_info.first_direction_unit;
      Tensor<1, dim> second_dir = parameters.ray_tracing_info.second_direction_unit;

      const double step_first_dir =
        parameters.ray_tracing_info.step_between_photons_first_direction;
      const double step_second_dir =
        parameters.ray_tracing_info.step_between_photons_second_direction;

      // Generation the photon insertion points.
      Point<dim> temp_point{};
      for (unsigned int n_first_dir = 0; n_first_dir < max_n_photon_first_dir;
           ++n_first_dir)
        {
          for (unsigned int n_second_dir = 0;
               n_second_dir < max_n_photon_second_dir;
               ++n_second_dir)
            {
              // Create the insertion point and emplace it.
              temp_point = [&]() {
                if constexpr (dim == 2)
                  return starting_insertion_point +
                         n_first_dir * step_first_dir * first_dir;

                if constexpr (dim == 3)
                  return starting_insertion_point +
                         n_first_dir * step_first_dir * first_dir +
                         n_second_dir * step_second_dir * second_dir;
              }();

              insertion_points_on_proc.push_back(temp_point);

              // Use the insertion point to create the property vector.
              std::vector<double> properties_of_one_photon(dim);

              properties_of_one_photon[0] = temp_point[0];
              properties_of_one_photon[1] = temp_point[1];
              if constexpr (dim == 3)
                properties_of_one_photon[2] = temp_point[2];

              photon_properties.push_back(properties_of_one_photon);
              properties_of_one_photon.clear();
            }
        }
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
  this->print_insertion_info(n_particles_to_insert_this_proc, pcout);
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
  for (const auto &p : points)
    {
      for (unsigned int d = 0; d < dim; ++d)
        oss << p[d] << " ";
      oss << "\n";
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
RayTracingSolver<dim>::finish_simulation()
{
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    this->computing_timer.print_summary();

  // Testing
  if (parameters.test.enabled)
    {
      // This needs to be coded later.
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
      for (auto cell_neighbor_list_iterator = cells_local_neighbor_list.begin();
           cell_neighbor_list_iterator != cells_local_neighbor_list.end();
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
              const Point<dim> photon_insertion_point = [&]() {
                if constexpr (dim == 2)
                  return Point<dim>(
                    {photon_properties[0], photon_properties[1]});

                if constexpr (dim == 3)
                  return Point<dim>({photon_properties[0],
                                     photon_properties[1],
                                     photon_properties[2]});
              }();

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
                  typename Particles::ParticleHandler<
                    dim>::particle_iterator_range particle_in_neighboring_cell =
                    particle_handler.particles_in_cell(
                      *current_neighboring_cell);

                  // This vector will be reused for every line-sphere
                  // intersection search. There is maximum 2 intersection points
                  // between a line and a sphere.
                  std::vector<Point<dim>> current_intersection_points;
                  current_intersection_points.reserve(2);

                  // Loop of the particles in the current neighboring cells
                  for (auto current_particle =
                         particle_in_neighboring_cell.begin();
                       current_particle != particle_in_neighboring_cell.end();
                       ++current_particle)
                    {
                      // Check if the current photon intersects the current
                      // particle
                      const Point<dim> particle_position =
                        current_particle->get_location();
                      const double particle_diameter =
                        current_particle->get_properties()
                          [DEMProperties::PropertiesIndex::dp];

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
                          const double distance_0 =
                            (photon_insertion_point -
                             current_intersection_points[0])
                              .norm();
                          const double distance_1 =
                            (photon_insertion_point -
                             current_intersection_points[1])
                              .norm();

                          if (distance_0 < distance_1)
                            {
                              new_distance = distance_0;
                              new_closest_point =
                                current_intersection_points[0];
                            }
                          else
                            {
                              new_distance = distance_1;
                              new_closest_point =
                                current_intersection_points[1];
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
                          photon_intersection_points_map[current_photon
                                                           ->get_id()] = {
                            new_distance, new_closest_point, current_photon};
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
                                new_distance,
                                new_closest_point,
                                current_photon};
                            }
                        }
                    }
                }
              // Even if we are removing photon at the end on this pseudo-time
              // step, we move the particle here since we are looping on each of
              // them in this loop. This way, we don't need to loop on each of
              // them at the end.
              const Point<dim> new_photon_location =
                current_photon_location + photon_displacement_vector;
              current_photon->set_location(new_photon_location);
            }
        }

      for (auto cell_neighbor_list_iterator = cells_ghost_neighbor_list.begin();
           cell_neighbor_list_iterator != cells_ghost_neighbor_list.end();
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

              Point<dim> photon_insertion_point;
              // Current photon insertion point.
              photon_insertion_point[0] = photon_properties[0];
              photon_insertion_point[1] = photon_properties[1];
              if constexpr (dim == 3)
                photon_insertion_point[2] = photon_properties[2];


              // Loop over the neighboring cells of the main cell. This includes
              // the main cell itselft.
              auto &cell_neighbor_list = *cell_neighbor_list_iterator;
              for (auto current_neighboring_cell = cell_neighbor_list.begin();
                   current_neighboring_cell != cell_neighbor_list.end();
                   ++current_neighboring_cell)
                {
                  // Particle iterators in the current neighboring cell
                  typename Particles::ParticleHandler<
                    dim>::particle_iterator_range particle_in_neighboring_cell =
                    particle_handler.particles_in_cell(
                      *current_neighboring_cell);

                  // This vector will be reused for every line-sphere
                  // intersection search. There is maximum 2 intersection points
                  // between a line and a sphere.
                  std::vector<Point<dim>> current_intersection_points;
                  current_intersection_points.reserve(2);

                  // Loop of the particles in the current neighboring cells
                  for (auto current_particle =
                         particle_in_neighboring_cell.begin();
                       current_particle != particle_in_neighboring_cell.end();
                       ++current_particle)
                    {
                      // Check if the current photon intersects the current
                      // particle
                      const Point<dim> particle_position =
                        current_particle->get_location();
                      const double particle_diameter =
                        current_particle->get_properties()
                          [DEMProperties::PropertiesIndex::dp];

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
                          const double distance_0 =
                            (photon_insertion_point -
                             current_intersection_points[0])
                              .norm();
                          const double distance_1 =
                            (photon_insertion_point -
                             current_intersection_points[1])
                              .norm();

                          if (distance_0 < distance_1)
                            {
                              new_distance = distance_0;
                              new_closest_point =
                                current_intersection_points[0];
                            }
                          else
                            {
                              new_distance = distance_1;
                              new_closest_point =
                                current_intersection_points[1];
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
                          photon_intersection_points_map[current_photon
                                                           ->get_id()] = {
                            new_distance, new_closest_point, current_photon};
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
                                new_distance,
                                new_closest_point,
                                current_photon};
                            }
                        }
                    }
                }
            }
        }
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
      photon_intersection_points_map.clear();

      action_manager->reset_triggers();
    }
  write_output_results(total_intersection_points,
                       parameters.simulation_control.output_folder,
                       parameters.simulation_control.output_name);

  finish_simulation();
}


template class RayTracingSolver<2>;
template class RayTracingSolver<3>;
