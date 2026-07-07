// SPDX-FileCopyrightText: Copyright (c) 2020-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/manifolds.h>
#include <core/solutions_output.h>

#include <dem/data_containers.h>
#include <dem/dem.h>
#include <dem/dem_post_processing.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/input_parameter_inspection.h>
#include <dem/lagrangian_post_processing.h>
#include <dem/multiphysics_integrator.h>
#include <dem/output_force_torque_calculation.h>
#include <dem/read_checkpoint.h>
#include <dem/read_mesh.h>
#include <dem/set_insertion_method.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/utilities.h>
#include <dem/velocity_verlet_integrator.h>
#include <dem/write_checkpoint.h>

#include <deal.II/base/table_handler.h>

#include <deal.II/grid/grid_out.h>

#include <sys/stat.h>

#include <numbers>
#include <ranges>
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
                    TimerOutput::never,
                    TimerOutput::wall_times)
  , contact_build_number(0)
  , background_dh(triangulation)
  , size_distribution_object_container(
      parameters.lagrangian_physical_properties.particle_type_number)
  , is_packed_insertion_method(false)
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

  // If the output directory does not exist, create it
  if (this_mpi_process == 0)
    {
      if (stat(output_dir_name.c_str(), &buffer) != 0)
        {
          create_output_folder(output_dir_name);
        }
    }

  // Get the pointer of the only instance of the action manager
  action_manager = DEMActionManager::get_action_manager();

  // Set the simulation control as transient DEM
  simulation_control = std::make_shared<SimulationControlTransientDEM>(
    parameters.simulation_control);

  // Set up the load balancing parameters and attach the correct functions to
  // the signals inside the triangulation
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
  grid_motion_object = std::make_shared<::GridMotionBase<dim, dim>>(
    parameters.grid_motion, simulation_control->get_time_step());

  // Set up the solid objects
  setup_solid_objects();

  // Check whether periodic boundaries are present. If at least one periodic
  // boundary is found, initialize the information for all periodic boundaries.
  if (!parameters.boundary_conditions.periodic_direction.empty())
    periodic_boundaries_object.set_periodic_boundaries_information(
      parameters.boundary_conditions.periodic_direction);

  // Assign gravity
  g = parameters.lagrangian_physical_properties.g;

  // If this is a restart simulation
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
  maximum_particle_diameter = 0;
  setup_distributions(parameters.lagrangian_physical_properties,
                      size_distribution_object_container,
                      maximum_particle_diameter,
                      this_mpi_process,
                      pcout);

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

  // Set the insertion object type before the restart because the restart only
  // rebuilds the member of the insertion object
  insertion_object =
    set_insertion_type<dim, PropertiesIndex>(size_distribution_object_container,
                                             triangulation,
                                             parameters,
                                             maximum_particle_diameter,
                                             is_packed_insertion_method);

  if (is_packed_insertion_method)
    {
      disable_position_integration = true;
      parameters.model_parameters.particle_particle_contact_force_model =
        ParticleParticleContactForceModel::shift;
      parameters.model_parameters.particle_wall_contact_force_method =
        ParticleWallContactForceModel::shift;

      if (parameters.model_parameters.sparse_particle_contacts)
        {
          sparse_contacts_object =
            AdaptiveSparseContacts<dim, PropertiesIndex>();
          DEMActionManager::get_action_manager()
            ->set_sparse_contacts_disabled();
        }
    }

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
  // the find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(triangulation) -
              maximum_particle_diameter * 0.5),
             (parameters.model_parameters.dynamic_contact_search_factor *
              (parameters.model_parameters.neighborhood_threshold - 1) *
              maximum_particle_diameter * 0.5));

  // Find the smallest cell size and use this as the floating mesh mapping
  // criterion. The edge case comes when the cells are completely square/cubic.
  // In that case, every side of a cell is 2^-0.5 or 3^-0.5 times the
  // cell_diameter. We want to refresh the mapping each time a solid-object
  // passes through a cell or there will be late contact detection. Thus, we use
  // this value.
  smallest_solid_object_mapping_criterion = [&] {
    if constexpr (dim == 2) // 2^-0.5 * D_c,min
      return 0.5 * std::numbers::sqrt2 *
             GridTools::minimal_cell_diameter(triangulation);
    if constexpr (dim == 3) // 3^-0.5 * D_c,min
      return std::numbers::inv_sqrt3 *
             GridTools::minimal_cell_diameter(triangulation);
  }();

  // Set up background dof
  setup_background_dofs();

  // Set up the periodic boundaries (if PBC enabled)
  periodic_boundaries_object.map_periodic_cells(
    triangulation, periodic_boundaries_cells_information);

  // Set the combined_periodic_offsets to contact managers and particles contact
  // forces for periodic contact detection (if PBC enabled)
  contact_manager.set_combined_periodic_offsets(
    periodic_boundaries_object.get_combined_periodic_offsets());
  particle_particle_contact_force_object->set_combined_periodic_offsets(
    periodic_boundaries_object.get_combined_periodic_offsets());

  // Set the periodic offsets of the periodic boundary pairs for other classes
  for (const auto &pb_id :
       periodic_boundaries_object.get_periodic_directions() | std::views::keys)
    {
      particle_particle_contact_force_object->set_periodic_offset(
        periodic_boundaries_object.get_periodic_offset_distance(pb_id), pb_id);
    }

  // Set up the local and ghost cells (if ASC enabled)
  sparse_contacts_object.update_local_and_ghost_cell_set(background_dh);
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::setup_background_dofs()
{
  FE_Q<dim> background_fe(1);
  background_dh.distribute_dofs(background_fe);

  // Periodic nodes must be mapped, otherwise the disabling of contacts will not
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

      // Loop over the periodic boundary conditions, keyed by the principal
      // periodic boundary id (id0).
      for (auto const &[id0, id1] :
           parameters.boundary_conditions.periodic_neighbor_id)
        {
          const unsigned int direction =
            parameters.boundary_conditions.periodic_direction.at(id0);

          // Default boundaries contain information for periodic boundary
          // conditions that indicate id0 and id1 are 0 as default value To
          // ensure these default values are not parsed, only make the
          // periodicity constraints if id0 and id1 are distinct
          if (id0 != id1)
            DoFTools::make_periodicity_constraints(
              background_dh, id0, id1, direction, background_constraints);
        }

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

  // If the load balancing uses the sparse_contact object to calculate the
  // weight, make sure that the sparse contact object has all mobility status
  // correctly identified and refreshed.
  if (parameters.model_parameters.load_balance_method ==
      ModelParameters<dim>::LoadBalanceMethod::dynamic_with_sparse_contacts)
    sparse_contacts_object.identify_mobility_status(
      background_dh,
      particle_handler,
      triangulation.n_active_cells(),
      mpi_communicator);

  // Prepare the particle handler for the adaptation of the triangulation to the
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
    (simulation_control->get_iteration_number() %
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
  if ((simulation_control->get_iteration_number() %
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

  if ((simulation_control->get_iteration_number() %
       parameters.insertion_info.insertion_frequency) == 1 ||
      simulation_control->get_iteration_number() == 1)
    {
      insertion_object->insert(particle_handler, triangulation, parameters);

      action_manager->particle_insertion_step();
    }

  if (is_packed_insertion_method)
    update_previous_position();
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
          contact_manager.get_particle_floating_mesh_potentially_in_contact(),
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
    {
      this->pcout << std::defaultfloat;
      this->computing_timer.print_summary();
      this->pcout << std::scientific;
    }

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
          default:
            break;
        }
    }

  // Outputting force and torques over the boundaries
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
  const unsigned int iter        = simulation_control->get_iteration_number();
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

      // Set up background dofs to initiate right sized boundary_id vector
      Vector<float> boundary_id(background_dh.n_dofs());

      // Attach the boundary_id to the data_out_faces object
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
      write_post_processing_results<dim>(
        triangulation,
        grid_pvdhandler,
        background_dh,
        particle_handler,
        parameters,
        simulation_control->get_current_time(),
        simulation_control->get_iteration_number(),
        mpi_communicator,
        sparse_contacts_object);
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::report_statistics()
{
  // Update statistics on the contact list
  double number_of_list_built_since_last_log =
    double(contact_build_number) - contact_list.total;
  contact_list.max =
    std::max(number_of_list_built_since_last_log, contact_list.max);
  contact_list.min =
    std::min(number_of_list_built_since_last_log, contact_list.min);
  contact_list.total += number_of_list_built_since_last_log;
  contact_list.average = contact_list.total /
                         (simulation_control->get_iteration_number()) *
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
      this->pcout << std::defaultfloat;
      this->computing_timer.print_summary();
      this->pcout << std::scientific;
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

      // Updating the moment of inertia container
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
DEMSolver<dim, PropertiesIndex>::update_previous_position()
{
  for (auto cell : triangulation.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      // Particles in the cell
      typename Particles::ParticleHandler<dim>::particle_iterator_range
        particles_in_cell = particle_handler.particles_in_cell(cell);

      // If the main cell is not empty
      if (particles_in_cell.empty())
        continue;

      for (auto particle_in_cell = particles_in_cell.begin();
           particle_in_cell != particles_in_cell.end();
           ++particle_in_cell)
        {
          auto particle_properties = particle_in_cell->get_properties();

          Point<dim> particle_previous_position =
            particle_in_cell->get_location();

          // Since we don't want to create a PropertiesIndex only for this
          // insertion method, we use the velocity to store the previous
          // location.
          particle_properties[PropertiesIndex::v_x] =
            particle_previous_position[0];
          particle_properties[PropertiesIndex::v_y] =
            particle_previous_position[1];
          if constexpr (dim == 3)
            particle_properties[PropertiesIndex::v_z] =
              particle_previous_position[2];
        }
    }
}

template <int dim, typename PropertiesIndex>
void
DEMSolver<dim, PropertiesIndex>::clamp_displacement()
{
  const double max_disp = maximum_particle_diameter;

  for (auto particle = particle_handler.begin();
       particle != particle_handler.end();
       ++particle)
    {
      auto particle_properties = particle->get_properties();

      Point<dim> previous_position;
      previous_position[0] = particle_properties[PropertiesIndex::v_x];
      previous_position[1] = particle_properties[PropertiesIndex::v_y];
      if constexpr (dim == 3)
        previous_position[2] = particle_properties[PropertiesIndex::v_z];

      const Tensor<1, dim> displacement_tensor =
        particle->get_location() - previous_position;

      const double disp_norm = displacement_tensor.norm();

      // No movement
      if (std::isnan(disp_norm))
        continue;

      const unsigned int particle_id = particle->get_local_index();
      if (disp_norm > max_disp)
        {
          // Clamp position to at most max_disp from the previous position
          particle->set_location(previous_position +
                                 (max_disp / disp_norm) * displacement_tensor);
          displacement[particle_id] += max_disp;
        }
      else
        {
          displacement[particle_id] += disp_norm;
        }
    }
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

  report_cell_size_to_particle_diameter_ratio(triangulation,
                                              maximum_particle_diameter,
                                              pcout,
                                              mpi_communicator);

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
          if (simulation_control->get_iteration_number() == 0)
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

      if (is_packed_insertion_method)
        {
          unsigned int number_of_pp_contact_on_proc =
            particle_particle_contact_force_object->get_number_of_contacts();

          unsigned int number_of_pw_contact_on_proc =
            particle_wall_contact_force_object->get_number_of_contacts();

          const unsigned int total_number_of_pp_contacts =
            Utilities::MPI::sum(number_of_pp_contact_on_proc, mpi_communicator);

          const unsigned int total_number_of_pw_contacts =
            Utilities::MPI::sum(number_of_pw_contact_on_proc, mpi_communicator);
          if (total_number_of_pp_contacts + total_number_of_pw_contacts == 0)
            {
              pcout << "No contact detected. Exiting simulation." << std::endl;
              write_output_results();
              break;
            }
          if (simulation_control->is_verbose_iteration())
            pcout << std::endl
                  << "Total number of p-p contacts: "
                  << total_number_of_pp_contacts << std::endl
                  << "Total number of p-w contacts: "
                  << total_number_of_pw_contacts << std::endl
                  << std::endl;

          particle_particle_contact_force_object->reset_number_of_contacts();
          particle_wall_contact_force_object->reset_number_of_contacts();

          clamp_displacement();
          update_previous_position();
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
              (simulation_control->get_iteration_number() %
                 parameters.forces_torques.output_frequency ==
               0) &&
              (parameters.forces_torques.force_torque_verbosity ==
               Parameters::Verbosity::verbose))
            {
              write_forces_torques_output_locally(
                forces_boundary_information[simulation_control
                                              ->get_iteration_number()],
                torques_boundary_information[simulation_control
                                               ->get_iteration_number()]);
            }
        }

      // Post-processing if needed
      post_process_results();

      // Write checkpoint if needed
      if (checkpoint_controller.is_checkpoint_time_step(
            simulation_control->get_iteration_number()))
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
