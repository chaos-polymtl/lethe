// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/manifolds.h>
#include <core/solutions_output.h>

#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/dem_post_processing.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/input_parameter_inspection.h>
#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/multiphysics_integrator.h>
#include <dem/output_force_torque_calculation.h>
#include <dem/read_checkpoint.h>
#include <dem/read_mesh.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/write_checkpoint.h>

#include <deal.II/base/table_handler.h>

#include <deal.II/grid/grid_out.h>

#include <sys/stat.h>

#include <numbers>
#include <sstream>
#include <utility>

template <int dim, typename PropertiesIndex>
DEMSolver<dim, PropertiesIndex>::DEMSolver(
  DEMSolverParameters<dim> dem_parameters)
  : mpi_communicator(MPI_COMM_WORLD)
  , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator))
  , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator))
  , pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
  , parameters(std::move(dem_parameters))
  , checkpoint_controller(parameters.restart)
  , triangulation(this->mpi_communicator)
  , mapping(1)
  , particle_handler(triangulation, mapping, PropertiesIndex::n_properties)
  , computing_timer(this->mpi_communicator,
                    this->pcout,
                    TimerOutput::summary,
                    TimerOutput::wall_times)
  , contact_build_number(0)
  , background_dh(triangulation)
  , size_distribution_object_container(
      parameters.lagrangian_physical_properties.particle_type_number)
{}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_parameters()
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

  // Change the behavior of the timer for situations when you don't want outputs
  if (parameters.timer.type == Parameters::Timer::Type::none)
    computing_timer.disable_output();

  // Set the simulation control as transient DEM
  simulation_control = std::make_shared<SimulationControlTransientDEM>(
    parameters.simulation_control);

  // Setup load balancing parameters and attach the correct functions to the
  // signals inside the triangulation
  load_balancing.set_parameters(parameters.model_parameters);
  load_balancing.copy_references(simulation_control,
                                 triangulation,
                                 particle_handler,
                                 sparse_contacts_object);
  load_balancing.connect_weight_signals();

  // Set the adaptive sparse contacts parameters
  if (parameters.model_parameters.sparse_particle_contacts)
    {
      sparse_contacts_object.set_parameters(
        parameters.model_parameters.granular_temperature_threshold,
        parameters.model_parameters.solid_fraction_threshold,
        parameters.model_parameters.advect_particles);
    }

  // Set the distribution type and initialize the neighborhood threshold
  setup_distribution_type();

  if (this_mpi_process == 0)
    input_parameter_inspection(parameters,
                               pcout,
                               size_distribution_object_container);

  // Set the grid motion type
  grid_motion_object =
    std::make_shared<GridMotion<dim, dim>>(parameters.grid_motion,
                                           simulation_control->get_time_step());

  // Set up the solid objects
  setup_solid_objects();

  // Check if there are periodic boundaries
  for (unsigned int i_bc = 0;
       i_bc < parameters.boundary_conditions.bc_types.size();
       ++i_bc)
    {
      if (parameters.boundary_conditions.bc_types[i_bc] ==
          Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
        {
          periodic_boundaries_object.set_periodic_boundaries_information(
            parameters.boundary_conditions.periodic_boundary_0,
            parameters.boundary_conditions.periodic_direction);
          break;
        }
    }

  // Assign gravity/acceleration
  g = parameters.lagrangian_physical_properties.g;

  // If this is a start simulation
  if (parameters.restart.restart)
    {
      action_manager->restart_simulation();

      if (parameters.model_parameters.load_balance_method !=
          Parameters::Lagrangian::ModelParameters<dim>::LoadBalanceMethod::none)
        action_manager->load_balance_step();
    }

  // Disable position integration
  disable_position_integration =
    parameters.model_parameters.disable_position_integration;
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_distribution_type()
{
  // Use namespace and alias to make the code more readable
  using namespace Parameters::Lagrangian;
  LagrangianPhysicalProperties &lpp = parameters.lagrangian_physical_properties;

  maximum_particle_diameter = 0;
  for (unsigned int particle_type = 0; particle_type < lpp.particle_type_number;
       particle_type++)
    {
      switch (lpp.distribution_type.at(particle_type))
        {
          case SizeDistributionType::uniform:
            size_distribution_object_container[particle_type] =
              std::make_shared<UniformDistribution>(
                lpp.particle_average_diameter.at(particle_type));
            break;
          case SizeDistributionType::normal:
            size_distribution_object_container[particle_type] =
              std::make_shared<NormalDistribution>(
                lpp.particle_average_diameter.at(particle_type),
                lpp.particle_size_std.at(particle_type),
                lpp.seed_for_distributions[particle_type] + this_mpi_process);
            break;
          case SizeDistributionType::custom:
            size_distribution_object_container[particle_type] =
              std::make_shared<CustomDistribution>(
                lpp.particle_custom_diameter.at(particle_type),
                lpp.particle_custom_probability.at(particle_type),
                lpp.seed_for_distributions[particle_type] + this_mpi_process);
            break;
        }

      maximum_particle_diameter = std::max(
        maximum_particle_diameter,
        size_distribution_object_container[particle_type]->find_max_diameter());
    }

  neighborhood_threshold_squared =
    std::pow(parameters.model_parameters.neighborhood_threshold *
               maximum_particle_diameter,
             2);
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_solid_objects()
{
  // Set up solid objects and carry them in vectors
  for (unsigned int i_solid = 0;
       i_solid < parameters.solid_objects->number_solid_surfaces;
       ++i_solid)
    {
      solid_surfaces.push_back(std::make_shared<SerialSolid<dim - 1, dim>>(
        this->parameters.solid_objects->solid_surfaces[i_solid], i_solid));
    }

  for (unsigned int i_solid = 0;
       i_solid < parameters.solid_objects->number_solid_volumes;
       ++i_solid)
    {
      solid_volumes.push_back(std::make_shared<SerialSolid<dim, dim>>(
        this->parameters.solid_objects->solid_volumes[i_solid], i_solid));
    }

  // Resize the mesh info containers
  solid_surfaces_mesh_info.resize(solid_surfaces.size());
  solid_volumes_mesh_info.resize(solid_volumes.size());

  // Simulation has solid objects and resize the container
  if ((solid_surfaces.size() + solid_volumes.size()) > 0)
    {
      action_manager->set_solid_objects_enabled();
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_functions_and_pointers()
{
  contact_detection_iteration_check_function =
    set_contact_search_iteration_function();

  // Set insertion object type before the restart because the restart only
  // rebuilds the member of the insertion object
  insertion_object = set_insertion_type();

  // Setting chosen contact force, insertion and integration methods
  integrator_object = set_integrator_type();
  particle_particle_contact_force_object =
    set_particle_particle_contact_force_model<dim, PropertiesIndex>(parameters);
  particle_wall_contact_force_object =
    set_particle_wall_contact_force_model<dim, PropertiesIndex>(parameters);
}

template <int dim, typename PropertiesIndex>
std::function<void()>
DEMSolver<dim, PropertiesIndex>::set_contact_search_iteration_function()
{
  using namespace Parameters::Lagrangian;
  typename ModelParameters<dim>::ContactDetectionMethod
    contact_detection_method =
      parameters.model_parameters.contact_detection_method;

  switch (contact_detection_method)
    {
      case ModelParameters<dim>::ContactDetectionMethod::constant:
        return [&] { check_contact_search_iteration_constant(); };
      case ModelParameters<dim>::ContactDetectionMethod::dynamic:
        return [&] { check_contact_search_iteration_dynamic(); };
      default:
        throw(std::runtime_error("Invalid contact detection method."));
    }
}

template <int dim, typename PropertiesIndex>
std::shared_ptr<Insertion<dim, PropertiesIndex>>
DEMSolver<dim, PropertiesIndex>::set_insertion_type()
{
  using namespace Parameters::Lagrangian;
  typename InsertionInfo<dim>::InsertionMethod insertion_method =
    parameters.insertion_info.insertion_method;

  switch (insertion_method)
    {
      case InsertionInfo<dim>::InsertionMethod::file:
        {
          return std::make_shared<InsertionFile<dim, PropertiesIndex>>(
            size_distribution_object_container, triangulation, parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::list:
        {
          return std::make_shared<InsertionList<dim, PropertiesIndex>>(
            size_distribution_object_container, triangulation, parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::plane:
        {
          return std::make_shared<InsertionPlane<dim, PropertiesIndex>>(
            size_distribution_object_container, triangulation, parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::volume:
        {
          return std::make_shared<InsertionVolume<dim, PropertiesIndex>>(
            size_distribution_object_container,
            triangulation,
            parameters,
            maximum_particle_diameter);
        }
      default:
        throw(std::runtime_error("Invalid insertion method."));
    }
}

template <int dim, typename PropertiesIndex>
std::shared_ptr<Integrator<dim, PropertiesIndex>>
DEMSolver<dim, PropertiesIndex>::set_integrator_type()
{
  using namespace Parameters::Lagrangian;
  typename ModelParameters<dim>::IntegrationMethod integration_method =
    parameters.model_parameters.integration_method;

  switch (integration_method)
    {
      case ModelParameters<dim>::IntegrationMethod::velocity_verlet:
        return std::make_shared<
          VelocityVerletIntegrator<dim, PropertiesIndex>>();
      case ModelParameters<dim>::IntegrationMethod::explicit_euler:
        return std::make_shared<
          ExplicitEulerIntegrator<dim, PropertiesIndex>>();
      default:
        throw(std::runtime_error("Invalid integration method."));
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_triangulation_dependent_parameters()
{
  // Find the smallest contact search frequency criterion between (smallest
  // cell size - largest particle radius) and (security factor * (blob
  // diameter - 1) * the largest particle radius). This value is used in
  // find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(triangulation) -
              maximum_particle_diameter * 0.5),
             (parameters.model_parameters.dynamic_contact_search_factor *
              (parameters.model_parameters.neighborhood_threshold - 1) *
              maximum_particle_diameter * 0.5));

  // Find the smallest cell size and use this as the floating mesh mapping
  // criterion. The edge case comes when the cell are completely square/cubic.
  // In that case, every side of a cell are 2^-0.5 or 3^-0.5 times the
  // cell_diameter. We want to refresh the mapping each time the solid-objet
  // pass through a cell or there will be late contact detection. Thus, we use
  // this value.
  smallest_solid_object_mapping_criterion = [&] {
    if constexpr (dim == 2) // 2^-0.5 * D_c,min
      return 0.5 * std::numbers::sqrt2 *
             GridTools::minimal_cell_diameter(triangulation);
    if constexpr (dim == 3) // 3^-0.5 * D_c,min
      return std::numbers::inv_sqrt3 *
             GridTools::minimal_cell_diameter(triangulation);
  }();

  // Setup background dof
  setup_background_dofs();

  // Set up the periodic boundaries (if PBC enabled)
  periodic_boundaries_object.map_periodic_cells(
    triangulation, periodic_boundaries_cells_information);

  // Set the periodic offset to contact managers and particles contact forces
  // for periodic contact detection (if PBC enabled)
  contact_manager.set_periodic_offset(
    periodic_boundaries_object.get_periodic_offset_distance());
  particle_particle_contact_force_object->set_periodic_offset(
    periodic_boundaries_object.get_periodic_offset_distance());

  // Set up the local and ghost cells (if ASC enabled)
  sparse_contacts_object.update_local_and_ghost_cell_set(background_dh);
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_background_dofs()
{
  FE_Q<dim> background_fe(1);
  background_dh.distribute_dofs(background_fe);

  // Periodic nodes must be mapped otherwise the disabling of contacts will not
  // propagate the mobility status to the periodic nodes during iterations.
  // Identification of periodic nodes is done with the background constraints.
  // Those constraints are not used for any matrix assembly in DEM solver, this
  // approach comes from CFD-DEM coupling where void fraction constraints are
  // used to achieve the periodic node mapping.
  if (action_manager->check_periodic_boundaries_enabled() &&
      action_manager->check_sparse_contacts_enabled())
    {
      IndexSet locally_relevant_dofs =
        DoFTools::extract_locally_relevant_dofs(background_dh);

      background_constraints.clear();
      background_constraints.reinit(background_dh.locally_owned_dofs(),
                                    locally_relevant_dofs);

      DoFTools::make_periodicity_constraints(
        background_dh,
        parameters.boundary_conditions.periodic_boundary_0,
        parameters.boundary_conditions.periodic_boundary_1,
        parameters.boundary_conditions.periodic_direction,
        background_constraints);

      background_constraints.close();

      sparse_contacts_object.map_periodic_nodes(background_constraints);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::load_balance()
{
  load_balancing.check_load_balance_iteration();

  if (!action_manager->check_load_balance())
    return;

  TimerOutput::Scope t(this->computing_timer, "Load balancing");
  // Prepare particle handler for the adaptation of the triangulation to the
  // load
  particle_handler.prepare_for_coarsening_and_refinement();

  pcout << "-->Repartitioning triangulation" << std::endl;
  triangulation.repartition();

  // Unpack the particle handler after the mesh has been repartitioned
  particle_handler.unpack_after_coarsening_and_refinement();

  // If PBC are enabled, update the periodic cells
  periodic_boundaries_object.map_periodic_cells(
    triangulation, periodic_boundaries_cells_information);

  // If ASC is enabled, update the local and ghost cell set
  sparse_contacts_object.update_local_and_ghost_cell_set(background_dh);

  // Update neighbors of cells after load balance
  contact_manager.update_cell_neighbors(triangulation,
                                        periodic_boundaries_cells_information);

  boundary_cell_object.build(
    triangulation,
    parameters.floating_walls,
    parameters.boundary_conditions.outlet_boundaries,
    parameters.mesh.check_for_diamond_cells,
    parameters.mesh.expand_particle_wall_contact_search,
    pcout);

  // Update the boundary information (if grid motion)
  boundary_cell_object.update_boundary_info_after_grid_motion(
    updated_boundary_points_and_normal_vectors);

  const auto average_minimum_maximum_cells =
    Utilities::MPI::min_max_avg(triangulation.n_active_cells(),
                                mpi_communicator);

  const auto average_minimum_maximum_particles =
    Utilities::MPI::min_max_avg(particle_handler.n_locally_owned_particles(),
                                mpi_communicator);

  pcout << "Load balance finished " << std::endl;
  pcout
    << "Average, minimum and maximum number of particles on the processors are "
    << average_minimum_maximum_particles.avg << " , "
    << average_minimum_maximum_particles.min << " and "
    << average_minimum_maximum_particles.max << std::endl;
  pcout << "Minimum and maximum number of cells owned by the processors are "
        << average_minimum_maximum_cells.min << " and "
        << average_minimum_maximum_cells.max << std::endl;

  setup_background_dofs();
}

template <int dim, typename PropertiesIndex>
inline void
DEMSolver<dim, PropertiesIndex>::check_contact_search_iteration_dynamic()
{
  const bool parallel_update =
    (simulation_control->get_step_number() %
     parameters.model_parameters.contact_detection_frequency) == 0;
  find_particle_contact_detection_step<dim, PropertiesIndex>(
    particle_handler,
    simulation_control->get_time_step(),
    smallest_contact_search_criterion,
    mpi_communicator,
    displacement,
    parallel_update);
}

template <int dim, typename PropertiesIndex>
inline void
DEMSolver<dim, PropertiesIndex>::check_contact_search_iteration_constant()
{
  if ((simulation_control->get_step_number() %
       parameters.model_parameters.contact_detection_frequency) == 0)
    action_manager->contact_detection_step();
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::insert_particles()
{
  // If the insertion frequency is set to 0, then no particles are going
  // to be inserted
  if (parameters.insertion_info.insertion_frequency == 0)
    return;

  if ((simulation_control->get_step_number() %
       parameters.insertion_info.insertion_frequency) == 1 ||
      simulation_control->get_step_number() == 1)
    {
      insertion_object->insert(particle_handler, triangulation, parameters);

      action_manager->particle_insertion_step();
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::particle_wall_contact_force()
{
  // Particle-wall contact force
  particle_wall_contact_force_object->calculate_particle_wall_contact(
    contact_manager.get_particle_wall_in_contact(),
    simulation_control->get_time_step(),
    contact_outcome);

  // Particle-floating wall contact force
  if (parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_contact_force_object->calculate_particle_wall_contact(
        contact_manager.get_particle_floating_wall_in_contact(),
        simulation_control->get_time_step(),
        contact_outcome);
    }

  // Particle-solid objects contact force
  if (action_manager->check_solid_objects_enabled()) // until refactor
    {
      particle_wall_contact_force_object
        ->calculate_particle_solid_object_contact(
          contact_manager.get_particle_floating_mesh_in_contact(),
          simulation_control->get_time_step(),
          solid_surfaces,
          contact_outcome);
    }

  particle_point_line_contact_force_object
    .calculate_particle_point_contact_force(
      &contact_manager.get_particle_points_in_contact(),
      parameters.lagrangian_physical_properties,
      force);

  if constexpr (dim == 3)
    {
      particle_point_line_contact_force_object
        .calculate_particle_line_contact_force(
          &contact_manager.get_particle_lines_in_contact(),
          parameters.lagrangian_physical_properties,
          force);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::move_solid_objects()
{
  if (!action_manager->check_solid_objects_enabled())
    return;

  // Move the solid triangulations, previous time must be used here
  // instead of current time.
  for (auto &solid_object : solid_surfaces)
    solid_object->move_solid_triangulation(
      simulation_control->get_time_step(),
      simulation_control->get_previous_time());

  for (auto &solid_object : solid_volumes)
    solid_object->move_solid_triangulation(
      simulation_control->get_time_step(),
      simulation_control->get_previous_time());
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::update_temperature_solid_objects()
{
  // Previous time must be used here instead of current time, which is the time
  // for which we are doing calculations. The solid object is moved and its
  // temperature is updated before calculating the contact outcomes for the
  // current time step, so its properties should correspond to the
  // previous time step.
  if constexpr (std::is_same_v<PropertiesIndex,
                               DEM::DEMMPProperties::PropertiesIndex>)
    {
      if (!action_manager->check_solid_objects_enabled())
        return;

      for (auto &solid_object : solid_surfaces)
        solid_object->update_solid_temperature(
          simulation_control->get_previous_time());

      for (auto &solid_object : solid_volumes)
        solid_object->update_solid_temperature(
          simulation_control->get_previous_time());
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::finish_simulation()
{
  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::end)
    this->computing_timer.print_summary();

  // Testing
  if (parameters.test.enabled)
    {
      switch (parameters.test.test_type)
        {
          case Parameters::Testing::TestType::particles:
            {
              visualization_object.print_xyz(particle_handler,
                                             mpi_communicator,
                                             pcout);
              break;
            }
          case Parameters::Testing::TestType::mobility_status:
            {
              // Get mobility status vector sorted by cell id
              Vector<float> mobility_status(triangulation.n_active_cells());
              sparse_contacts_object.get_mobility_status_vector(
                mobility_status);

              // Output mobility status vector
              visualization_object.print_intermediate_format(mobility_status,
                                                             background_dh,
                                                             mpi_communicator);
              break;
            }
          case Parameters::Testing::TestType::subdomain:
            {
              // Get mobility status vector sorted by cell id
              Vector<float> subdomain(triangulation.n_active_cells());
              for (unsigned int i = 0; i < subdomain.size(); ++i)
                subdomain(i) = triangulation.locally_owned_subdomain();

              // Output subdomain vector
              visualization_object.print_intermediate_format(subdomain,
                                                             background_dh,
                                                             mpi_communicator);
              break;
            }
        }
    }

  // Outputting force and torques over boundary
  if (parameters.forces_torques.calculate_force_torque)
    {
      write_forces_torques_output_results(
        parameters.forces_torques.force_torque_output_name,
        parameters.forces_torques.output_frequency,
        triangulation.get_boundary_ids(),
        simulation_control->get_time_step(),
        forces_boundary_information,
        torques_boundary_information);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::write_output_results()
{
  TimerOutput::Scope t(this->computing_timer, "Output VTU");

  const std::string folder = parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    parameters.simulation_control.output_name;
  const unsigned int iter        = simulation_control->get_step_number();
  const double       time        = simulation_control->get_current_time();
  const unsigned int group_files = parameters.simulation_control.group_files;

  // Write particles
  Visualization<dim, PropertiesIndex> particle_data_out;
  particle_data_out.build_patches(particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler,
                            particle_data_out,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            mpi_communicator);

  if (simulation_control->get_output_boundaries())
    {
      DataOutFaces<dim> data_out_faces;

      // Setup background dofs to initiate right sized boundary_id vector
      Vector<float> boundary_id(background_dh.n_dofs());

      // Attach the boundary_id to data_out_faces object
      BoundaryPostprocessor<dim> boundary;
      data_out_faces.attach_dof_handler(background_dh);
      data_out_faces.add_data_vector(boundary_id, boundary);
      data_out_faces.build_patches();

      write_boundaries_vtu<dim>(
        data_out_faces, folder, time, iter, this->mpi_communicator);
    }
  if (parameters.post_processing.force_chains)
    {
      // Force chains visualization
      particles_force_chains_object =
        set_force_chains_contact_force_model<dim, PropertiesIndex>(parameters);
      particles_force_chains_object->write_force_chains(
        parameters,
        particles_pvdhandler_force_chains,
        this->mpi_communicator,
        folder,
        iter,
        time,
        contact_manager.get_local_adjacent_particles(),
        contact_manager.get_ghost_adjacent_particles());
    }

  // Write all solid objects
  for (const auto &solid_object : solid_surfaces)
    solid_object->write_output_results(simulation_control);

  for (const auto &solid_object : solid_volumes)
    solid_object->write_output_results(simulation_control);
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::post_process_results()
{
  if (parameters.post_processing.lagrangian_post_processing_enabled &&
      simulation_control->is_output_iteration())
    {
      write_post_processing_results<dim>(triangulation,
                                         grid_pvdhandler,
                                         background_dh,
                                         particle_handler,
                                         parameters,
                                         simulation_control->get_current_time(),
                                         simulation_control->get_step_number(),
                                         mpi_communicator,
                                         sparse_contacts_object);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::report_statistics()
{
  // Update statistics on contact list
  double number_of_list_built_since_last_log =
    double(contact_build_number) - contact_list.total;
  contact_list.max =
    std::max(number_of_list_built_since_last_log, contact_list.max);
  contact_list.min =
    std::min(number_of_list_built_since_last_log, contact_list.min);
  contact_list.total += number_of_list_built_since_last_log;
  contact_list.average = contact_list.total /
                         (simulation_control->get_step_number()) *
                         simulation_control->get_log_frequency();

  // Calculate statistics on the particles
  statistics translational_kinetic_energy = calculate_granular_statistics<
    dim,
    PropertiesIndex,
    DEM::dem_statistic_variable::translational_kinetic_energy>(
    particle_handler, mpi_communicator);
  statistics rotational_kinetic_energy = calculate_granular_statistics<
    dim,
    PropertiesIndex,
    DEM::dem_statistic_variable::rotational_kinetic_energy>(particle_handler,
                                                            mpi_communicator);
  statistics velocity =
    calculate_granular_statistics<dim,
                                  PropertiesIndex,
                                  DEM::dem_statistic_variable::velocity>(
      particle_handler, mpi_communicator);
  statistics omega =
    calculate_granular_statistics<dim,
                                  PropertiesIndex,
                                  DEM::dem_statistic_variable::omega>(
      particle_handler, mpi_communicator);

  if (this_mpi_process == 0)
    {
      TableHandler report;

      report.declare_column("Variable");
      report.declare_column("Min");
      report.declare_column("Max");
      report.declare_column("Average");
      report.declare_column("Total");
      add_statistics_to_table_handler("Contact list generation",
                                      contact_list,
                                      report);
      add_statistics_to_table_handler("Velocity magnitude", velocity, report);
      add_statistics_to_table_handler("Angular velocity magnitude",
                                      omega,
                                      report);
      add_statistics_to_table_handler("Translational kinetic energy",
                                      translational_kinetic_energy,
                                      report);
      add_statistics_to_table_handler("Rotational kinetic energy",
                                      rotational_kinetic_energy,
                                      report);



      report.set_scientific("Min", true);
      report.set_scientific("Max", true);
      report.set_scientific("Average", true);
      report.set_scientific("Total", true);

      report.write_text(std::cout, dealii::TableHandler::org_mode_table);
    }

  // Timer output
  if (parameters.timer.type == Parameters::Timer::Type::iteration)
    {
      this->computing_timer.print_summary();
    }
}

template <int dim, typename PropertiesIndex>
inline void
DEMSolver<dim, PropertiesIndex>::sort_particles_into_subdomains_and_cells()
{
  particle_handler.sort_particles_into_subdomains_and_cells();

  // Exchange ghost particles
  particle_handler.exchange_ghost_particles(true);

  // Resize the displacement, force and torque containers only if the particles
  // have changed subdomains
  if (action_manager->check_resize_containers())
    {
      unsigned int number_of_particles =
        particle_handler.get_max_local_particle_index();
      // Resize displacement container
      displacement.resize(number_of_particles);
      // Resize outcome containers
      contact_outcome.resize_interaction_containers(number_of_particles);
      MOI.resize(number_of_particles);

      // Updating moment of inertia container
      for (auto &particle : particle_handler)
        {
          auto particle_properties = particle.get_properties();
          MOI[particle.get_local_index()] =
            0.1 * particle_properties[PropertiesIndex::mass] *
            particle_properties[PropertiesIndex::dp] *
            particle_properties[PropertiesIndex::dp];
        }
    }

  // Always reset the displacement values since we are doing a search detection
  std::ranges::fill(displacement, 0.);
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::solve()
{
  // Set up the parameters
  setup_parameters();

  // Reading mesh
  read_mesh(parameters.mesh,
            action_manager->check_restart_simulation(),
            pcout,
            triangulation,
            parameters.boundary_conditions);

  // Set up functions and pointers according to parameters
  setup_functions_and_pointers();

  // Read checkpoint if needed
  read_checkpoint(computing_timer,
                  parameters,
                  simulation_control,
                  particles_pvdhandler,
                  grid_pvdhandler,
                  triangulation,
                  particle_handler,
                  insertion_object,
                  solid_surfaces,
                  checkpoint_controller);

  // Set up the various parameters that need the triangulation
  setup_triangulation_dependent_parameters();

  // Build the mapping of the cell neighbors
  contact_manager.execute_cell_neighbors_search(
    triangulation, periodic_boundaries_cells_information);

  // Find boundary cells with faces
  boundary_cell_object.build(
    triangulation,
    parameters.floating_walls,
    parameters.boundary_conditions.outlet_boundaries,
    parameters.mesh.check_for_diamond_cells,
    parameters.mesh.expand_particle_wall_contact_search,
    pcout);

  // DEM engine iterator
  while (simulation_control->integrate())
    {
      simulation_control->print_progression(pcout);
      if (simulation_control->is_verbose_iteration())
        report_statistics();

      // Move grid and update the boundary information (if grid motion)
      grid_motion_object->move_grid(triangulation);
      boundary_cell_object.update_boundary_info_after_grid_motion(
        updated_boundary_points_and_normal_vectors);

      // Insert particle if needed
      insert_particles();

      // Load balancing (if load balancing enabled and if needed)
      load_balance();

      // Check for contact search according to the contact detection method
      contact_detection_iteration_check_function();

      // Check if solid object needs to be mapped with background mesh
      // (if solid object)
      find_floating_mesh_mapping_step(smallest_solid_object_mapping_criterion,
                                      this->solid_surfaces);

      // Map solid objects if the action was triggered (if solid object)
      if (action_manager->check_solid_object_search())
        {
          // Store information about floating mesh/background mesh
          // intersection
          for (unsigned int i_solid = 0; i_solid < solid_surfaces.size();
               ++i_solid)
            {
              solid_surfaces_mesh_info[i_solid] =
                solid_surfaces[i_solid]->map_solid_in_background_triangulation(
                  triangulation);
            }

          for (unsigned int i_solid = 0; i_solid < solid_volumes.size();
               ++i_solid)
            {
              solid_volumes_mesh_info[i_solid] =
                solid_volumes[i_solid]->map_solid_in_background_triangulation(
                  triangulation);
            }
        }

      // Execute contact search if the action was triggered
      if (action_manager->check_contact_search())
        {
          // Particles displacement if passing through a periodic boundary
          // (if PBC enabled)
          periodic_boundaries_object.execute_particles_displacement(
            particle_handler, periodic_boundaries_cells_information);

          sort_particles_into_subdomains_and_cells();

          // Compute cell mobility (if ASC enabled)
          sparse_contacts_object.identify_mobility_status(
            background_dh,
            particle_handler,
            triangulation.n_active_cells(),
            mpi_communicator);

          // Execute broad search by filling containers of particle-particle
          // contact pair candidates and containers of particle-wall
          // contact pair candidates
          contact_manager.execute_particle_particle_broad_search(
            particle_handler, sparse_contacts_object);

          contact_manager.execute_particle_wall_broad_search(
            particle_handler,
            boundary_cell_object,
            solid_surfaces_mesh_info,
            parameters.floating_walls,
            simulation_control->get_current_time(),
            sparse_contacts_object);

          // Update contacts, remove replicates and add new contact pairs
          // to the contact containers when particles are exchanged between
          // processors
          contact_manager.update_contacts();

          // Updates the iterators to particles in local-local contact
          // containers
          contact_manager.update_local_particles_in_cells(particle_handler);

          // Execute fine search by updating particle-particle contact
          // containers according to the neighborhood threshold
          contact_manager.execute_particle_particle_fine_search(
            neighborhood_threshold_squared);

          // Execute fine search by updating particle-wall contact
          // containers according to the neighborhood threshold
          contact_manager.execute_particle_wall_fine_search(
            parameters.floating_walls,
            simulation_control->get_current_time(),
            neighborhood_threshold_squared);

          // Updating number of contact builds
          contact_build_number++;
        }
      else
        {
          particle_handler.update_ghost_particles();
        }

      // Particle-particle contact force
      particle_particle_contact_force_object
        ->calculate_particle_particle_contact(
          contact_manager.get_local_adjacent_particles(),
          contact_manager.get_ghost_adjacent_particles(),
          contact_manager.get_local_local_periodic_adjacent_particles(),
          contact_manager.get_local_ghost_periodic_adjacent_particles(),
          contact_manager.get_ghost_local_periodic_adjacent_particles(),
          simulation_control->get_time_step(),
          contact_outcome);

      // Update the boundary points and vectors (if grid motion)
      // We have to update the positions of the points on boundary faces and
      // their normal vectors here. The update_contacts deletes the
      // particle-wall contact candidate if it exists in the contact list.
      // As a result, when we update the points on boundary faces and their
      // normal vectors, update_contacts deletes it from the output of broad
      // search, and they are not updated in the contact force calculations.
      grid_motion_object
        ->update_boundary_points_and_normal_vectors_in_contact_list(
          contact_manager.get_particle_wall_in_contact(),
          updated_boundary_points_and_normal_vectors);

      // Move solid objects (if solid object)
      move_solid_objects();

      // Update solid objects temperatures
      update_temperature_solid_objects();

      // Particle-wall contact force
      particle_wall_contact_force();

      // Integration of temperature for multiphysic DEM
      if constexpr (std::is_same_v<PropertiesIndex,
                                   DEM::DEMMPProperties::PropertiesIndex>)
        {
          integrate_temperature<dim, PropertiesIndex>(
            particle_handler,
            simulation_control->get_time_step(),
            contact_outcome.heat_transfer_rate,
            std::vector<double>(force.size()));
        }

      // Integration of force and velocity for new location of particles
      // The half step is calculated at the first iteration

      if (!disable_position_integration)
        {
          if (simulation_control->get_step_number() == 0)
            {
              integrator_object->integrate_half_step_location(
                particle_handler,
                g,
                simulation_control->get_time_step(),
                torque,
                force,
                MOI);
            }
          else
            {
              integrator_object->integrate(particle_handler,
                                           g,
                                           simulation_control->get_time_step(),
                                           torque,
                                           force,
                                           MOI,
                                           triangulation,
                                           sparse_contacts_object);
            }
        }

      // Visualization
      if (simulation_control->is_output_iteration())
        write_output_results();

      // Log the contact statistics if the parameter is enabled
      if (parameters.post_processing.particle_wall_collision_statistics)
        {
          log_collision_data<dim, PropertiesIndex>(
            parameters,
            contact_manager.get_particle_wall_in_contact(),
            simulation_control->get_current_time(),
            ongoing_collision_log,
            collision_event_log);
        }

      // Calculation of forces and torques if needed
      if (parameters.forces_torques.calculate_force_torque)
        {
          if ((this_mpi_process == 0) &&
              (simulation_control->get_step_number() %
                 parameters.forces_torques.output_frequency ==
               0) &&
              (parameters.forces_torques.force_torque_verbosity ==
               Parameters::Verbosity::verbose))
            {
              write_forces_torques_output_locally(
                forces_boundary_information[simulation_control
                                              ->get_step_number()],
                torques_boundary_information[simulation_control
                                               ->get_step_number()]);
            }
        }

      // Post-processing if needed
      post_process_results();

      // Write checkpoint if needed
      if (checkpoint_controller.is_checkpoint_time_step(
            simulation_control->get_step_number()))
        {
          write_checkpoint(computing_timer,
                           parameters,
                           simulation_control,
                           particles_pvdhandler,
                           grid_pvdhandler,
                           triangulation,
                           particle_handler,
                           insertion_object,
                           solid_surfaces,
                           pcout,
                           mpi_communicator,
                           checkpoint_controller);
        }

      // Reset all trigger flags
      action_manager->reset_triggers();
    }

  // Write particle-wall collision statistics file if enabled
  if (parameters.post_processing.particle_wall_collision_statistics)
    write_collision_stats(parameters, collision_event_log, mpi_communicator);

  finish_simulation();
}

template class DEMSolver<2, DEM::DEMProperties::PropertiesIndex>;
template class DEMSolver<3, DEM::DEMProperties::PropertiesIndex>;
template class DEMSolver<2, DEM::DEMMPProperties::PropertiesIndex>;
template class DEMSolver<3, DEM::DEMMPProperties::PropertiesIndex>;
