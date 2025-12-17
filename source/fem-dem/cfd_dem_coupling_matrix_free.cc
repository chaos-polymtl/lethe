// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>
#include <core/solutions_output.h>

#include <solvers/postprocessing_cfd.h>
#include <solvers/postprocessing_velocities.h>

#include <dem/dem_post_processing.h>
#include <dem/explicit_euler_integrator.h>
#include <dem/find_contact_detection_step.h>
#include <dem/insertion.h>
#include <dem/insertion_file.h>
#include <dem/insertion_list.h>
#include <dem/insertion_plane.h>
#include <dem/insertion_volume.h>
#include <dem/particle_handler_conversion.h>
#include <dem/set_particle_particle_contact_force_model.h>
#include <dem/set_particle_wall_contact_force_model.h>
#include <dem/velocity_verlet_integrator.h>
#include <fem-dem/cfd_dem_coupling_matrix_free.h>
#include <fem-dem/fluid_dynamics_vans_matrix_free_operators.h>
#include <fem-dem/postprocessing_cfd_dem.h>


// Constructor for the class CFD-DEM class
template <int dim>
CFDDEMMatrixFree<dim>::CFDDEMMatrixFree(CFDDEMSimulationParameters<dim> &param)
  : FluidDynamicsVANSMatrixFree<dim>(param)
{}

template <int dim>
void
CFDDEMMatrixFree<dim>::setup_distribution_type()
{
  // Use namespace and alias to make the code more readable
  using namespace Parameters::Lagrangian;
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(this->mpi_communicator);
  LagrangianPhysicalProperties &lpp =
    dem_parameters.lagrangian_physical_properties;

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
                lpp.seed_for_distributions[particle_type] + this_mpi_process,
                lpp.diameter_min_cutoff.at(particle_type),
                lpp.diameter_max_cutoff.at(particle_type));
            break;
          case SizeDistributionType::lognormal:
            size_distribution_object_container[particle_type] =
              std::make_shared<LogNormalDistribution>(
                lpp.particle_average_diameter.at(particle_type),
                lpp.particle_size_std.at(particle_type),
                lpp.seed_for_distributions[particle_type] + this_mpi_process,
                lpp.diameter_min_cutoff.at(particle_type),
                lpp.diameter_max_cutoff.at(particle_type));
            break;
          case SizeDistributionType::custom:
            size_distribution_object_container[particle_type] =
              std::make_shared<CustomDistribution>(
                lpp.particle_custom_diameter.at(particle_type),
                lpp.particle_custom_probability.at(particle_type),
                lpp.seed_for_distributions[particle_type] + this_mpi_process);
            break;
        }
      size_distribution_object_container[particle_type]
        ->print_psd_declaration_string(particle_type, this->pcout);

      maximum_particle_diameter = std::max(
        maximum_particle_diameter,
        size_distribution_object_container[particle_type]->find_max_diameter());
    }

  neighborhood_threshold_squared =
    std::pow(dem_parameters.model_parameters.neighborhood_threshold *
               maximum_particle_diameter,
             2);
}

template <int dim>
void
CFDDEMMatrixFree<dim>::dem_setup_parameters()
{
  coupling_frequency =
    this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency;

  dem_action_manager = DEMActionManager::get_action_manager();

  // Initialize DEM Parameters
  dem_parameters.lagrangian_physical_properties =
    this->cfd_dem_simulation_parameters.dem_parameters
      .lagrangian_physical_properties;
  g = dem_parameters.lagrangian_physical_properties.g;
  dem_parameters.boundary_conditions =
    this->cfd_dem_simulation_parameters.dem_parameters.boundary_conditions;
  dem_parameters.insertion_info =
    this->cfd_dem_simulation_parameters.dem_parameters.insertion_info;
  dem_parameters.floating_walls =
    this->cfd_dem_simulation_parameters.dem_parameters.floating_walls;
  dem_parameters.model_parameters =
    this->cfd_dem_simulation_parameters.dem_parameters.model_parameters;
  dem_parameters.simulation_control =
    this->cfd_dem_simulation_parameters.dem_parameters.simulation_control;
  dem_parameters.post_processing =
    this->cfd_dem_simulation_parameters.dem_parameters.post_processing;
  dem_parameters.mesh = this->cfd_dem_simulation_parameters.dem_parameters.mesh;
  dem_parameters.restart =
    this->cfd_dem_simulation_parameters.dem_parameters.restart;
  size_distribution_object_container.resize(
    dem_parameters.lagrangian_physical_properties.particle_type_number);

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  // Setup load balancing parameters
  load_balancing.set_parameters(dem_parameters.model_parameters);
  load_balancing.copy_references(this->simulation_control,
                                 *parallel_triangulation,
                                 this->particle_handler,
                                 sparse_contacts_object);

  // Attach the correct functions to the signals inside the triangulation
  load_balancing.connect_weight_signals();

  if (dem_parameters.model_parameters.sparse_particle_contacts)
    {
      dem_action_manager->set_sparse_contacts_enabled();
      sparse_contacts_object.set_parameters(
        dem_parameters.model_parameters.granular_temperature_threshold,
        dem_parameters.model_parameters.solid_fraction_threshold,
        dem_parameters.model_parameters.advect_particles);
    }

  setup_distribution_type();

  // Calculate the Rayleigh critical time step ratio
  rayleigh_time_step = DBL_MAX;

  for (unsigned int i = 0;
       i < dem_parameters.lagrangian_physical_properties.particle_type_number;
       ++i)
    {
      double youngs_modulus = dem_parameters.lagrangian_physical_properties
                                .youngs_modulus_particle[i];
      double poisson_ratio =
        dem_parameters.lagrangian_physical_properties.poisson_ratio_particle[i];
      double density =
        dem_parameters.lagrangian_physical_properties.density_particle[i];

      double shear_modulus = youngs_modulus / (2.0 * (1.0 + poisson_ratio));

      double min_diameter =
        size_distribution_object_container.at(i)->find_min_diameter();

      rayleigh_time_step =
        std::min(M_PI_2 * min_diameter * sqrt(density / shear_modulus) /
                   (0.1631 * poisson_ratio + 0.8766),
                 rayleigh_time_step);
    }

  update_dem_time_step();

  // Check if there are periodic boundaries
  for (unsigned int i_bc = 0;
       i_bc < dem_parameters.boundary_conditions.bc_types.size();
       ++i_bc)
    {
      if (dem_parameters.boundary_conditions.bc_types[i_bc] ==
          Parameters::Lagrangian::BCDEM::BoundaryType::periodic)
        {
          dem_action_manager->set_periodic_boundaries_enabled();

          periodic_boundaries_object.set_periodic_boundaries_information(
            dem_parameters.boundary_conditions.periodic_boundary_0,
            dem_parameters.boundary_conditions.periodic_direction);

          break;
        }
    }

  insertion_object = set_insertion_type();
  // Initialize the total contact list counter
  integrator_object = set_integrator_type();
  particle_particle_contact_force_object =
    set_particle_particle_contact_force_model<
      dim,
      DEM::CFDDEMProperties::PropertiesIndex>(
      this->cfd_dem_simulation_parameters.dem_parameters);

  // Initialize the contact search counter
  contact_search_total_number = 0;
}


template <int dim>
std::shared_ptr<Insertion<dim, DEM::CFDDEMProperties::PropertiesIndex>>
CFDDEMMatrixFree<dim>::set_insertion_type()
{
  using namespace Parameters::Lagrangian;
  typename InsertionInfo<dim>::InsertionMethod insertion_method =
    dem_parameters.insertion_info.insertion_method;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  switch (insertion_method)
    {
      case InsertionInfo<dim>::InsertionMethod::file:
        {
          return std::make_shared<
            InsertionFile<dim, DEM::CFDDEMProperties::PropertiesIndex>>(
            size_distribution_object_container,
            *parallel_triangulation,
            dem_parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::list:
        {
          return std::make_shared<
            InsertionList<dim, DEM::CFDDEMProperties::PropertiesIndex>>(
            size_distribution_object_container,
            *parallel_triangulation,
            dem_parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::plane:
        {
          return std::make_shared<
            InsertionPlane<dim, DEM::CFDDEMProperties::PropertiesIndex>>(
            size_distribution_object_container,
            *parallel_triangulation,
            dem_parameters);
        }
      case InsertionInfo<dim>::InsertionMethod::volume:
        {
          return std::make_shared<
            InsertionVolume<dim, DEM::CFDDEMProperties::PropertiesIndex>>(
            size_distribution_object_container,
            *parallel_triangulation,
            dem_parameters,
            maximum_particle_diameter);
        }
      default:
        throw(std::runtime_error("Invalid insertion method."));
    }
}

template <int dim>
std::shared_ptr<Integrator<dim, DEM::CFDDEMProperties::PropertiesIndex>>
CFDDEMMatrixFree<dim>::set_integrator_type()
{
  using namespace Parameters::Lagrangian;
  typename ModelParameters<dim>::IntegrationMethod integration_method =
    dem_parameters.model_parameters.integration_method;

  switch (integration_method)
    {
      case ModelParameters<dim>::IntegrationMethod::velocity_verlet:
        return std::make_shared<
          VelocityVerletIntegrator<dim,
                                   DEM::CFDDEMProperties::PropertiesIndex>>();
      case ModelParameters<dim>::IntegrationMethod::explicit_euler:
        return std::make_shared<
          ExplicitEulerIntegrator<dim,
                                  DEM::CFDDEMProperties::PropertiesIndex>>();
      default:
        throw(std::runtime_error("Invalid integration method."));
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::initialize_dem_parameters()
{
  this->pcout << "Initializing DEM parameters" << std::endl;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);

  // Set up the local and ghost cells (if ASC enabled)
  sparse_contacts_object.update_local_and_ghost_cell_set(
    this->particle_projector.dof_handler);

  particle_wall_contact_force_object = set_particle_wall_contact_force_model<
    dim,
    DEM::CFDDEMProperties::PropertiesIndex>(
    this->cfd_dem_simulation_parameters.dem_parameters);

  // Finding the smallest contact search frequency criterion between (smallest
  // cell size - largest particle radius) and (security factor * (blob diameter
  // - 1) *  the largest particle radius). This value is used in
  // find_contact_detection_frequency function
  smallest_contact_search_criterion =
    std::min((GridTools::minimal_cell_diameter(*this->triangulation) -
              maximum_particle_diameter * 0.5),
             (dem_parameters.model_parameters.dynamic_contact_search_factor *
              (dem_parameters.model_parameters.neighborhood_threshold - 1) *
              maximum_particle_diameter * 0.5));

  // Remap periodic cells (if PBC enabled)
  periodic_boundaries_object.map_periodic_cells(
    *parallel_triangulation, periodic_boundaries_cells_information);

  // Set the periodic offset to contact managers and particles contact forces
  // for periodic contact detection (if PBC enabled)
  contact_manager.set_periodic_offset(
    periodic_boundaries_object.get_periodic_offset_distance());
  particle_particle_contact_force_object->set_periodic_offset(
    periodic_boundaries_object.get_periodic_offset_distance());

  // Find cell neighbors
  contact_manager.execute_cell_neighbors_search(
    *parallel_triangulation, periodic_boundaries_cells_information);

  // Find boundary cells with faces
  boundary_cell_object.build(
    *parallel_triangulation,
    dem_parameters.floating_walls,
    dem_parameters.boundary_conditions.outlet_boundaries,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .check_for_diamond_cells,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh
      .expand_particle_wall_contact_search,
    this->pcout);

  sort_particles_into_subdomains_and_cells();

  // Remap periodic nodes after setup of dofs (If ASC and PBC)
  if (dem_action_manager->check_periodic_boundaries_enabled() &&
      dem_action_manager->check_sparse_contacts_enabled())
    {
      sparse_contacts_object.map_periodic_nodes(
        this->particle_projector.void_fraction_constraints);
    }

  this->pcout << "Finished initializing DEM parameters" << std::endl
              << "DEM time-step is " << dem_time_step << " s" << std::endl;
}

template <int dim>
void
CFDDEMMatrixFree<dim>::read_dem()
{
  this->pcout << "Reading DEM checkpoint" << std::endl;

  std::string prefix =
    this->cfd_dem_simulation_parameters.void_fraction->dem_file_name;

  // Load checkpoint controller
  std::string checkpoint_controller_object_filename =
    prefix + ".checkpoint_controller";
  std::ifstream iss_checkpoint_controller_obj(
    checkpoint_controller_object_filename);
  boost::archive::text_iarchive ia_checkpoint_controller_obj(
    iss_checkpoint_controller_obj, boost::archive::no_header);

  unsigned int checkpoint_id;
  ia_checkpoint_controller_obj >> checkpoint_id;

  // New prefix for the remaining files
  prefix = prefix + "_" + Utilities::int_to_string(checkpoint_id);

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  // Create a temporary particle_handler with DEM properties
  Particles::ParticleHandler<dim> temporary_particle_handler(
    *this->triangulation,
    this->particle_mapping,
    DEM::DEMProperties::n_properties);

  ia >> temporary_particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      try
        {
          parallel_triangulation->load(filename.c_str());

          // Deserialize particles have the triangulation has been read
          temporary_particle_handler.deserialize();
        }
      catch (...)
        {
          AssertThrow(false,
                      ExcMessage("Cannot open snapshot mesh file or read the"
                                 "triangulation stored there."));
        }

      // Fill the existing particle handler using the temporary one
      // This is done during the dynamic cast for the convert_particle_handler
      // function which requires a parallel::distributed::triangulation
      convert_particle_handler<dim,
                               DEM::DEMProperties::PropertiesIndex,
                               DEM::CFDDEMProperties::PropertiesIndex>(
        *parallel_triangulation,
        temporary_particle_handler,
        this->particle_handler);
    }
  else
    {
      throw std::runtime_error(
        "VANS equations currently do not support triangulations other than parallel::distributed");
    }

  this->pcout << "Finished reading DEM checkpoint" << std::endl
              << this->particle_handler.n_global_particles()
              << " particles are in the simulation" << std::endl;

  write_dem_output_results();
}

template <int dim>
std::vector<OutputStructTableHandler>
CFDDEMMatrixFree<dim>::gather_tables()
{
  std::vector<OutputStructTableHandler> table_output_structs;

  const Parameters::PostProcessing post_processing =
    this->simulation_parameters.post_processing;
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder;
  std::string suffix = ".checkpoint";
  if (post_processing.calculate_phase_volumes)
    table_output_structs.emplace_back(
      this->table_phase_volumes,
      prefix + post_processing.phase_volumes_output_name + suffix);
  return table_output_structs;
}

template <int dim>
void
CFDDEMMatrixFree<dim>::write_checkpoint()
{
  using VectorType = LinearAlgebra::distributed::Vector<double>;

  TimerOutput::Scope timer(this->computing_timer, "Write checkpoint");

  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  std::string prefix_grid      = prefix + "_lagrangian_postprocessing";
  if (Utilities::MPI::this_mpi_process(this->mpi_communicator) == 0)
    {
      this->simulation_control->save(prefix);
      this->pvdhandler.save(prefix);
      particles_pvdhandler.save(prefix_particles);

      if (this->dem_parameters.post_processing
            .lagrangian_post_processing_enabled)
        grid_pvdhandler.save(prefix_grid);

      if (this->simulation_parameters.flow_control.enable_flow_control)
        this->flow_control.save(prefix);
    }

  std::ostringstream            oss;
  boost::archive::text_oarchive oa(oss, boost::archive::no_header);
  oa << this->particle_handler;

  // Write additional particle information for deserialization
  std::string   particle_filename = prefix + ".particles";
  std::ofstream output(particle_filename.c_str());
  output << oss.str() << std::endl;

  std::vector<const VectorType *> sol_set_transfer;
  sol_set_transfer.push_back(&(*this->present_solution));
  for (unsigned int i = 0; i < this->previous_solutions->size(); ++i)
    {
      sol_set_transfer.push_back(&(*this->previous_solutions)[i]);
    }

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<const VectorType *> av_set_transfer =
        this->average_velocities->save(prefix);

      // Insert average velocities vectors into the set transfer vector
      sol_set_transfer.insert(sol_set_transfer.end(),
                              av_set_transfer.begin(),
                              av_set_transfer.end());
    }

  // Prepare for Serialization
  SolutionTransfer<dim, VectorType> system_trans_vectors(*this->dof_handler);
  system_trans_vectors.prepare_for_serialization(sol_set_transfer);

  // Prepare particle handler for serialization
  this->particle_handler.prepare_for_serialization();

  // Void Fraction
  std::vector<const VectorType *> vf_set_transfer;
  vf_set_transfer.push_back(&this->particle_projector.void_fraction_solution);
  for (unsigned int i = 0;
       i < this->particle_projector.void_fraction_previous_solution.size();
       ++i)
    {
      vf_set_transfer.push_back(
        &this->particle_projector.void_fraction_previous_solution[i]);
    }

  this->multiphysics->write_checkpoint();

  // Prepare for Serialization
  SolutionTransfer<dim, VectorType> vf_system_trans_vectors(
    this->particle_projector.dof_handler);
  vf_system_trans_vectors.prepare_for_serialization(vf_set_transfer);

  if (auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation))
    {
      std::string triangulationName = prefix + ".triangulation";
      parallel_triangulation->save(prefix + ".triangulation");
    }
  // Serialize all post-processing tables that are currently used
  // Serialize the post-processing tables that are additional in this solver
  const std::vector<OutputStructTableHandler> &table_output_structs_add =
    this->gather_tables();
  serialize_tables_vector(table_output_structs_add, this->mpi_communicator);
  // Serialize the default post-processing tables that are members of
  // NavierStokesBase
  const std::vector<OutputStructTableHandler> &table_output_structs =
    NavierStokesBase<dim, VectorType, IndexSet>::gather_tables();
  serialize_tables_vector(table_output_structs, this->mpi_communicator);
}

template <int dim>
void
CFDDEMMatrixFree<dim>::read_checkpoint()
{
  TimerOutput::Scope timer(this->computing_timer, "Read checkpoint");
  using VectorType = LinearAlgebra::distributed::Vector<double>;
  std::string prefix =
    this->simulation_parameters.simulation_control.output_folder +
    this->simulation_parameters.restart_parameters.filename;
  std::string prefix_particles = prefix + "_particles";
  std::string prefix_grid      = prefix + "_lagrangian_postprocessing";

  this->simulation_control->read(prefix);
  this->pvdhandler.read(prefix);
  particles_pvdhandler.read(prefix_particles);

  if (this->dem_parameters.post_processing.lagrangian_post_processing_enabled)
    {
      grid_pvdhandler.read(prefix_grid);
    }

  // Gather particle serialization information
  std::string   particle_filename = prefix + ".particles";
  std::ifstream input(particle_filename.c_str());
  AssertThrow(input, ExcFileNotOpen(particle_filename));

  std::string buffer;
  std::getline(input, buffer);
  std::istringstream            iss(buffer);
  boost::archive::text_iarchive ia(iss, boost::archive::no_header);

  ia >> this->particle_handler;

  const std::string filename = prefix + ".triangulation";
  std::ifstream     in(filename.c_str());
  if (!in)
    AssertThrow(false,
                ExcMessage(
                  std::string(
                    "You are trying to restart a previous computation, "
                    "but the restart file <") +
                  filename + "> does not appear to exist!"));

  try
    {
      if (auto parallel_triangulation =
            dynamic_cast<parallel::distributed::Triangulation<dim> *>(
              this->triangulation.get()))
        parallel_triangulation->load(filename.c_str());
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("Cannot open snapshot mesh file or read the "
                             "triangulation stored there."));
    }

  this->setup_dofs();

  // TODO BB
  // Remap periodic nodes after setup of dofs
  if (dem_action_manager->check_periodic_boundaries_enabled() &&
      dem_action_manager->check_sparse_contacts_enabled())
    {
      sparse_contacts_object.map_periodic_nodes(
        this->particle_projector.void_fraction_constraints);
    }

  // Velocity Vectors
  std::vector<VectorType *> x_system(1 + this->previous_solutions->size());

  VectorType distributed_system(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  x_system[0] = &(distributed_system);

  std::vector<VectorType> distributed_previous_solutions;

  distributed_previous_solutions.reserve(this->previous_solutions->size());

  for (unsigned int i = 0; i < this->previous_solutions->size(); ++i)
    {
      distributed_previous_solutions.emplace_back(
        VectorType(this->locally_owned_dofs,
                   this->locally_relevant_dofs,
                   this->mpi_communicator));
      x_system[i + 1] = &distributed_previous_solutions[i];
    }

  SolutionTransfer<dim, VectorType> system_trans_vectors(*this->dof_handler);

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      std::vector<VectorType *> sum_vectors =
        this->average_velocities->read(prefix);

      x_system.insert(x_system.end(), sum_vectors.begin(), sum_vectors.end());
    }

  system_trans_vectors.deserialize(x_system);

  *this->present_solution = distributed_system;
  for (unsigned int i = 0; i < this->previous_solutions->size(); ++i)
    {
      (*this->previous_solutions)[i] = distributed_previous_solutions[i];
    }

  // Void Fraction Vectors
  std::vector<VectorType *> vf_system(
    1 + this->particle_projector.previous_void_fraction.size());

  VectorType vf_distributed_system(
    this->particle_projector.locally_owned_dofs,
    this->particle_projector.locally_relevant_dofs,
    this->mpi_communicator);

  vf_system[0] = &(vf_distributed_system);

  std::vector<VectorType> vf_distributed_previous_solutions;

  vf_distributed_previous_solutions.reserve(
    this->particle_projector.previous_void_fraction.size());

  for (unsigned int i = 0;
       i < this->particle_projector.previous_void_fraction.size();
       ++i)
    {
      vf_distributed_previous_solutions.emplace_back(
        VectorType(this->particle_projector.locally_owned_dofs,
                   this->particle_projector.locally_relevant_dofs,
                   this->mpi_communicator));
      vf_system[i + 1] = &vf_distributed_previous_solutions[i];
    }

  SolutionTransfer<dim, VectorType> vf_system_trans_vectors(
    this->particle_projector.dof_handler);

  vf_system_trans_vectors.deserialize(vf_system);

  this->particle_projector.void_fraction_solution = vf_distributed_system;

#ifndef LETHE_USE_LDV
  // We also wish the Trilinos solution to be updated.
  convert_vector_dealii_to_trilinos(
    this->particle_projector.void_fraction_locally_relevant,
    this->particle_projector.void_fraction_solution);
#endif

  for (unsigned int i = 0;
       i < this->particle_projector.previous_void_fraction.size();
       ++i)
    {
      this->particle_projector.void_fraction_previous_solution[i] =
        vf_distributed_previous_solutions[i];

#ifndef LETHE_USE_LDV
      // We also wish the Trilinos solution to be updated.
      convert_vector_dealii_to_trilinos(
        this->particle_projector.previous_void_fraction[i],
        this->particle_projector.void_fraction_previous_solution[i]);
#endif
    }

  if (this->simulation_parameters.flow_control.enable_flow_control)
    {
      this->flow_control.read(prefix);
    }
  this->multiphysics->read_checkpoint();

  // Deserialize particles have the triangulation has been read
  this->particle_handler.deserialize();

  // Deserialize all post-processing tables that are currently used
  // Deserialize the post-processing tables that are particular to this solver
  std::vector<OutputStructTableHandler> table_output_structs_add =
    this->gather_tables();
  deserialize_tables_vector(table_output_structs_add, this->mpi_communicator);

  // Deserialize the default post-processing tables that are members of
  // NavierStokesBase
  std::vector<OutputStructTableHandler> table_output_structs =
    NavierStokesBase<dim, VectorType, IndexSet>::gather_tables();
  deserialize_tables_vector(table_output_structs, this->mpi_communicator);
}

template <int dim>
void
CFDDEMMatrixFree<dim>::check_contact_detection_method(unsigned int counter)
{
  // Use namespace and alias to make the code more readable
  using namespace Parameters::Lagrangian;
  Parameters::Lagrangian::ModelParameters<dim> &model_parameters =
    dem_parameters.model_parameters;

  switch (model_parameters.contact_detection_method)
    {
      case ModelParameters<dim>::ContactDetectionMethod::constant:
        {
          if ((counter % model_parameters.contact_detection_frequency) == 0)
            dem_action_manager->contact_detection_step();
          break;
        }
      case ModelParameters<dim>::ContactDetectionMethod::dynamic:
        {
          double dt =
            this->simulation_control->get_time_step() /
            this->cfd_dem_simulation_parameters.cfd_dem.coupling_frequency;

          find_particle_contact_detection_step<
            dim,
            DEM::CFDDEMProperties::PropertiesIndex>(
            this->particle_handler,
            dt,
            smallest_contact_search_criterion,
            this->mpi_communicator,
            displacement);
          break;
        }
      default:
        break;
    }
}



template <int dim>
void
CFDDEMMatrixFree<dim>::load_balance()
{
  load_balancing.check_load_balance_iteration();

  // If not a load balance iteration, exit the function
  if (!dem_action_manager->check_load_balance())
    return;

  // Otherwise, we do not support load balancing at the present time, throw and
  // exit gracefully.
  AssertThrow(false, ExcMessage("Load balancing is currently not supported"));
}

template <int dim>
void
CFDDEMMatrixFree<dim>::add_fluid_particle_interaction()
{
  /// Reference to torque vector from contact outcomes
  std::vector<Tensor<1, 3>> &torque = contact_outcome.torque;

  /// Reference to force vector from contact outcomes
  std::vector<Tensor<1, 3>> &force = contact_outcome.force;

  for (auto particle = this->particle_handler.begin();
       particle != this->particle_handler.end();
       ++particle)
    {
      auto particle_properties = particle->get_properties();

      types::particle_index particle_id = particle->get_local_index();

      force[particle_id][0] +=
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_two_way_coupling_x] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_one_way_coupling_x] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::fem_drag_x];
      force[particle_id][1] +=
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_two_way_coupling_y] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_one_way_coupling_y] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::fem_drag_y];
      force[particle_id][2] +=
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_two_way_coupling_z] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::
                              fem_force_one_way_coupling_z] +
        particle_properties[DEM::CFDDEMProperties::PropertiesIndex::fem_drag_z];
      torque[particle_id][0] += particle_properties
        [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_x];
      torque[particle_id][1] += particle_properties
        [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_y];
      torque[particle_id][2] += particle_properties
        [DEM::CFDDEMProperties::PropertiesIndex::fem_torque_z];
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::insert_particles()
{
  // If the insertion frequency is set to 0, then no particles are going
  // to be inserted in the CFD-DEM simulation and the function returns
  if (dem_parameters.insertion_info.insertion_frequency == 0)
    return;

  const auto parallel_triangulation =
    dynamic_cast<parallel::distributed::Triangulation<dim> *>(
      &*this->triangulation);
  if ((this->simulation_control->get_step_number() %
       dem_parameters.insertion_info.insertion_frequency) == 1 ||
      this->simulation_control->get_step_number() == 1)
    {
      insertion_object->insert(this->particle_handler,
                               *parallel_triangulation,
                               dem_parameters);

      dem_action_manager->particle_insertion_step();

      // Sort particles after insertion
      sort_particles_into_subdomains_and_cells();
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::particle_wall_contact_force()
{
  // Particle-wall contact force
  particle_wall_contact_force_object->calculate_particle_wall_contact(
    contact_manager.get_particle_wall_in_contact(),
    dem_time_step,
    contact_outcome);

  // Particle-floating wall contact force
  if (dem_parameters.floating_walls.floating_walls_number > 0)
    {
      particle_wall_contact_force_object->calculate_particle_wall_contact(
        contact_manager.get_particle_floating_wall_in_contact(),
        dem_time_step,
        contact_outcome);
    }

  particle_point_line_contact_force_object
    .calculate_particle_point_contact_force(
      &contact_manager.get_particle_points_in_contact(),
      dem_parameters.lagrangian_physical_properties,
      contact_outcome.force);

  if constexpr (dim == 3)
    {
      particle_point_line_contact_force_object
        .calculate_particle_line_contact_force(
          &contact_manager.get_particle_lines_in_contact(),
          dem_parameters.lagrangian_physical_properties,
          contact_outcome.force);
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::write_dem_output_results()
{
  const std::string folder = dem_parameters.simulation_control.output_folder;
  const std::string particles_solution_name =
    dem_parameters.simulation_control.output_name + "_particles";
  const unsigned int iter = this->simulation_control->get_step_number();
  const double       time = this->simulation_control->get_current_time();
  const unsigned int group_files =
    dem_parameters.simulation_control.group_files;

  // Write particles
  Visualization<dim, DEM::CFDDEMProperties::PropertiesIndex> particle_data_out;
  particle_data_out.build_patches(this->particle_handler,
                                  properties_class.get_properties_name());

  write_vtu_and_pvd<0, dim>(particles_pvdhandler,
                            particle_data_out,
                            folder,
                            particles_solution_name,
                            time,
                            iter,
                            group_files,
                            this->mpi_communicator);
}


template <int dim>
void
CFDDEMMatrixFree<dim>::report_particle_statistics()
{
  const auto this_mpi_process =
    Utilities::MPI::this_mpi_process(this->mpi_communicator);
  // Update statistics on contact list
  double number_of_list_built_since_last_log =
    double(contact_search_total_number) - contact_list.total;
  contact_list.max =
    std::max(number_of_list_built_since_last_log, contact_list.max);
  contact_list.min =
    std::min(number_of_list_built_since_last_log, contact_list.min);
  contact_list.total   = double(contact_search_total_number);
  contact_list.average = contact_list.total /
                         (this->simulation_control->get_step_number()) *
                         this->simulation_control->get_log_frequency();

  // Calculate statistics on the particles
  statistics translational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::CFDDEMProperties::PropertiesIndex,
    DEM::dem_statistic_variable::translational_kinetic_energy>(
    this->particle_handler, this->mpi_communicator);
  statistics rotational_kinetic_energy = calculate_granular_statistics<
    dim,
    DEM::CFDDEMProperties::PropertiesIndex,
    DEM::dem_statistic_variable::rotational_kinetic_energy>(
    this->particle_handler, this->mpi_communicator);
  statistics velocity =
    calculate_granular_statistics<dim,
                                  DEM::CFDDEMProperties::PropertiesIndex,
                                  DEM::dem_statistic_variable::velocity>(
      this->particle_handler, this->mpi_communicator);
  statistics omega =
    calculate_granular_statistics<dim,
                                  DEM::CFDDEMProperties::PropertiesIndex,
                                  DEM::dem_statistic_variable::omega>(
      this->particle_handler, this->mpi_communicator);

  if (this_mpi_process == 0)
    {
      TableHandler report;

      std::vector<std::string> column_names{
        "Variable", "Min", "Max", "Average", "Total"};

      for (const std::string &column_name : column_names)
        report.declare_column(column_name);

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

      // Only for Min, Max, Average and Total columns
      for (unsigned int i = 1; i < column_names.size(); ++i)
        {
          report.set_scientific(column_names[i], true);
          report.set_precision(column_names[i],
                               this->simulation_control->get_log_precision());
        }

      announce_string(this->pcout, "Particle statistics");
      report.write_text(std::cout, dealii::TableHandler::org_mode_table);
    }

  if (dem_parameters.post_processing.lagrangian_post_processing_enabled &&
      this->simulation_control->is_output_iteration())
    {
      const auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation);

      write_post_processing_results<dim>(
        *parallel_triangulation,
        grid_pvdhandler,
        *this->dof_handler,
        this->particle_handler,
        dem_parameters,
        this->simulation_control->get_current_time(),
        this->simulation_control->get_step_number(),
        this->mpi_communicator,
        sparse_contacts_object);
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::postprocess_fd(bool first_iteration)
{
  this->FluidDynamicsMatrixFree<dim>::postprocess_fd(first_iteration);

  // Visualization
  if (this->simulation_control->is_output_iteration())
    {
      write_dem_output_results();
    }
  if (!first_iteration)
    {
      this->pcout
        << "---------------------------------------------------------------"
        << std::endl;
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::postprocess_cfd_dem()
{
  // Calculate total volume of fluid and solid
  if (this->simulation_parameters.post_processing.calculate_phase_volumes)
    {
      TimerOutput::Scope t(this->computing_timer, "total_volume_calculation");
      double             total_volume_fluid, total_volume_particles;
      std::tie(total_volume_fluid, total_volume_particles) =
        calculate_fluid_and_particle_volumes(
          this->particle_projector.dof_handler,
          this->particle_projector.void_fraction_locally_relevant,
          *this->cell_quadrature,
          *this->mapping);
      this->table_phase_volumes.add_value(
        "time", this->simulation_control->get_current_time());
      this->table_phase_volumes.add_value("total-volume-fluid",
                                          total_volume_fluid);
      this->table_phase_volumes.add_value("total-volume-particles",
                                          total_volume_particles);
      if (this->simulation_parameters.post_processing.verbosity ==
          Parameters::Verbosity::verbose)
        {
          this->pcout << "Total volume of fluid: "
                      << std::setprecision(
                           this->simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           total_volume_fluid
                      << " m^3" << std::endl;
          this->pcout << "Total volume of particles: "
                      << std::setprecision(
                           this->simulation_control->get_log_precision())
                      << this->simulation_parameters.physical_properties_manager
                             .get_density_scale() *
                           total_volume_particles
                      << " m^3" << std::endl;
        }

      // Output pressure drop to a text file from processor 0
      if ((this->simulation_control->get_step_number() %
             this->simulation_parameters.post_processing.output_frequency ==
           0) &&
          this->this_mpi_process == 0)
        {
          std::string filename =
            this->simulation_parameters.simulation_control.output_folder +
            this->simulation_parameters.post_processing
              .phase_volumes_output_name +
            ".dat";
          std::ofstream output(filename.c_str());
          table_phase_volumes.set_precision("time", 12);
          table_phase_volumes.set_precision("total-volume-fluid", 12);
          table_phase_volumes.set_precision("total-volume-particles", 12);
          this->table_phase_volumes.write_text(output);
        }
    }
}


template <int dim>
void
CFDDEMMatrixFree<dim>::dynamic_flow_control()
{
  if (this->simulation_parameters.flow_control.enable_flow_control &&
      this->simulation_parameters.simulation_control.method !=
        Parameters::SimulationControl::TimeSteppingMethod::steady)
    {
      // Calculate the average velocity according to the void fraction
      unsigned int flow_direction =
        this->simulation_parameters.flow_control.flow_direction;
      double average_velocity = calculate_average_velocity(
        *this->dof_handler,
        this->particle_projector.dof_handler,
        *this->present_solution,
        this->particle_projector.void_fraction_solution,
        flow_direction,
        *this->cell_quadrature,
        *this->mapping);

      // Calculate the beta force for fluid
      this->flow_control.calculate_beta(
        average_velocity,
        this->simulation_control->get_time_step(),
        this->simulation_control->get_step_number());

      Tensor<1, 3> beta_particle;
      if (this->simulation_parameters.flow_control.enable_beta_particle)
        {
          // Calculate the beta for particles and add it to force tensor g
          AssertThrow(
            dem_parameters.lagrangian_physical_properties
                .particle_type_number == 1,
            ExcMessage(std::string(
              "Beta force calculation for particles is not implemented for multiple particle types.")));

          double fluid_density =
            this->simulation_parameters.physical_properties_manager
              .get_density_scale();
          double particle_density =
            dem_parameters.lagrangian_physical_properties.density_particle[0];
          beta_particle =
            this->flow_control.get_beta_particles(fluid_density,
                                                  particle_density);
          g = dem_parameters.lagrangian_physical_properties.g + beta_particle;
        }

      // Showing results
      if (this->simulation_parameters.flow_control.verbosity ==
            Parameters::Verbosity::verbose &&
          this->simulation_control->get_step_number() > 0 &&
          this->this_mpi_process == 0)
        {
          announce_string(this->pcout, "Flow control summary");
          this->pcout << "Fluid space-average velocity: " << average_velocity
                      << std::endl;
          this->pcout << "Fluid beta force: "
                      << this->flow_control.get_beta()[flow_direction]
                      << std::endl;
          this->pcout << "Particle beta force: "
                      << beta_particle[flow_direction] << std::endl;
        }
    }
}


template <int dim>
inline void
CFDDEMMatrixFree<dim>::sort_particles_into_subdomains_and_cells()
{
  this->particle_handler.sort_particles_into_subdomains_and_cells();

  // Exchange ghost particles
  this->particle_handler.exchange_ghost_particles(true);

  // Resize the displacement, force and torque containers only if the particles
  // have changed subdomains
  if (dem_action_manager->check_resize_containers())
    {
      unsigned int number_of_particles =
        this->particle_handler.get_max_local_particle_index();
      // Resize displacement container
      displacement.resize(number_of_particles);
      // Resize outcome containers
      contact_outcome.resize_interaction_containers(number_of_particles);
      MOI.resize(number_of_particles);

      // Updating moment of inertia container
      for (auto &particle : this->particle_handler)
        {
          auto particle_properties = particle.get_properties();
          MOI[particle.get_local_index()] =
            0.1 *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::mass] *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp] *
            particle_properties[DEM::CFDDEMProperties::PropertiesIndex::dp];
        }
    }

  // Always reset the displacement values since we are doing a search detection
  std::ranges::fill(displacement, 0.);

  this->particle_handler.exchange_ghost_particles(true);
}

template <int dim>
void
CFDDEMMatrixFree<dim>::dem_iterator(unsigned int counter)
{
  // dem_contact_build carries out the particle-particle and particle-wall
  // broad and fine searches, sort_particles_into_subdomains_and_cells, and
  // exchange_ghost
  dem_contact_build(counter);

  // Particle-particle contact force
  particle_particle_contact_force_object->calculate_particle_particle_contact(
    contact_manager.get_local_adjacent_particles(),
    contact_manager.get_ghost_adjacent_particles(),
    contact_manager.get_local_local_periodic_adjacent_particles(),
    contact_manager.get_local_ghost_periodic_adjacent_particles(),
    contact_manager.get_ghost_local_periodic_adjacent_particles(),
    dem_time_step,
    contact_outcome);

  // Particles-walls contact force:
  particle_wall_contact_force();

  // Add fluid-particle interaction force to the force container
  add_fluid_particle_interaction();

  // The cell average velocities and accelerations are updated
  // from the fully computed forces at step 0, but we do not use the
  // mobility status to disabled contacts at this step. The reason
  // is that the current mobility status are computed for the last
  // DEM time step (from the last CFD time step) are the force of
  // the fluid may have significantly changed the particulate
  // agitation.
  // Update the cell average velocities and accelerations
  sparse_contacts_object.update_average_velocities_acceleration(
    this->particle_handler, g, contact_outcome.force, dem_time_step);

  // Integration correction step (after force calculation)
  // In the first step, we have to obtain location of particles at half-step
  // time
  // TODO do all DEM time step at first CFD time step are half step?
  if (this->simulation_control->get_step_number() == 0)
    {
      integrator_object->integrate_half_step_location(this->particle_handler,
                                                      g,
                                                      dem_time_step,
                                                      contact_outcome.torque,
                                                      contact_outcome.force,
                                                      MOI);
    }
  else
    {
      const auto parallel_triangulation =
        dynamic_cast<parallel::distributed::Triangulation<dim> *>(
          &*this->triangulation);
      integrator_object->integrate(this->particle_handler,
                                   g,
                                   dem_time_step,
                                   contact_outcome.torque,
                                   contact_outcome.force,
                                   MOI,
                                   *parallel_triangulation,
                                   sparse_contacts_object);
    }

  auto dem_current_time =
    (this->simulation_control->get_current_time()) + (dem_time_step * counter);

  // Log the contact statistics if the parameter is enabled
  if (dem_parameters.post_processing.particle_wall_collision_statistics)
    {
      log_collision_data<dim, DEM::CFDDEMProperties::PropertiesIndex>(
        dem_parameters,
        contact_manager.get_particle_wall_in_contact(),
        dem_current_time,
        ongoing_collision_log,
        collision_event_log);
    }

  dem_action_manager->reset_triggers();
}

template <int dim>
void
CFDDEMMatrixFree<dim>::dem_contact_build(unsigned int counter)
{
  // If this is not the last DEM iteration before the next CFD iteration, check
  // if a contact detection step is necessary
  if (counter != (coupling_frequency - 1))
    check_contact_detection_method(counter);

  // Otherwise, force a contact search at the last DEM iteration before a CFD
  // iteration to ensure that the particles are adequately located before
  // calculating the coupling between particle and fluid.
  else
    dem_action_manager->last_dem_of_cfddem_iteration_step();

  // Sort particles in cells
  if (dem_action_manager->check_contact_search())
    {
      this->pcout << "DEM contact search at dem step " << counter << std::endl;
      contact_search_counter++;

      // Execute periodic boundaries (if PBC enabled)
      periodic_boundaries_object.execute_particles_displacement(
        this->particle_handler, periodic_boundaries_cells_information);

      // Sort particles into subdomains and cells and reset the vectors that
      // the local particle if is their index
      sort_particles_into_subdomains_and_cells();

      // Identify the mobility status of particles (if ASC enabled)
      sparse_contacts_object.identify_mobility_status(
        this->particle_projector.dof_handler,
        this->particle_handler,
        (*this->triangulation).n_active_cells(),
        this->mpi_communicator);

      // Execute broad search by filling containers of particle-particle
      // contact pair candidates and containers of particle-wall
      // contact pair candidates
      contact_manager.execute_particle_particle_broad_search(
        this->particle_handler, sparse_contacts_object);

      contact_manager.execute_particle_wall_broad_search(
        this->particle_handler,
        boundary_cell_object,
        solid_surfaces_mesh_info,
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        sparse_contacts_object);

      // Update contacts, remove replicates and add new contact pairs
      // to the contact containers when particles are exchanged between
      // processors
      contact_manager.update_contacts();

      // Updates the iterators to particles in local-local contact
      // containers
      contact_manager.update_local_particles_in_cells(this->particle_handler);

      // Execute fine search by updating particle-particle contact
      // containers regards the neighborhood threshold
      contact_manager.execute_particle_particle_fine_search(
        neighborhood_threshold_squared);

      // Execute fine search by updating particle-wall contact containers
      // regards the neighborhood threshold
      contact_manager.execute_particle_wall_fine_search(
        dem_parameters.floating_walls,
        this->simulation_control->get_current_time(),
        neighborhood_threshold_squared);
    }
  else
    {
      this->particle_handler.update_ghost_particles();
    }
}

template <int dim>
void
CFDDEMMatrixFree<dim>::solve()
{
  this->computing_timer.enter_subsection("Read mesh, manifolds and particles");

  read_mesh_and_manifolds(
    *this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  dem_setup_parameters();

  // Reading DEM start file information
  if (this->cfd_dem_simulation_parameters.void_fraction->read_dem == true &&
      !this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
         .restart)
    read_dem();

  this->computing_timer.leave_subsection("Read mesh, manifolds and particles");


  this->setup_dofs();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  if (this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart)
    dem_action_manager->restart_simulation();

  // Initialize the DEM parameters and generate the ghost particles
  initialize_dem_parameters();

  // Calculate first instance of void fraction once particles are set up
  if (!dem_action_manager->check_restart_simulation())
    this->particle_projector.initialize_void_fraction(
      this->simulation_control->get_current_time());

  // Output the solution after initializing the void fraction
  if (!this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
         .restart)
    {
      this->postprocess_fd(true);
      this->multiphysics->postprocess(true);
      if (this->simulation_control->is_output_iteration())
        this->write_output_results(*this->present_solution);
    }

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);

      // We allow the physics to update their boundary conditions
      // according to their own parameters
      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      // Insert particle if needed
      insert_particles();

      this->dynamic_flow_control();

      if (!this->simulation_control->is_at_start())
        {
          this->refine_mesh();
        }

      // We calculate the void fraction and the particle-fluid interaction using
      // the particle projector.
      this->particle_projector.calculate_void_fraction(
        this->simulation_control->get_current_time());

      if (time_stepping_is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->computing_timer.enter_subsection("Calculate time derivatives");

          this->calculate_time_derivative_previous_solutions();
          this->time_derivative_previous_solutions.update_ghost_values();
          this->system_operator->evaluate_time_derivative_previous_solutions(
            this->time_derivative_previous_solutions);
          this->evaluate_time_derivative_void_fraction();
          this->computing_timer.leave_subsection("Calculate time derivatives");

          if (this->simulation_parameters.flow_control.enable_flow_control)
            this->system_operator->update_beta_force(
              this->flow_control.get_beta());
        }


      this->particle_projector.calculate_particle_fluid_forces_projection(
        this->cfd_dem_simulation_parameters.cfd_dem,
        *this->dof_handler,
        *this->present_solution,
        *this->previous_solutions,
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties.g,
        NavierStokesScratchData<dim>(
          this->simulation_control,
          this->simulation_parameters.physical_properties_manager,
          *this->fe,
          *this->cell_quadrature,
          *this->mapping,
          *this->face_quadrature));

      // The base matrix-free operator is not aware of the various VANS
      // coupling terms. We must do a cast here to ensure that the operator is
      // of the right type.
      if (auto mf_operator = dynamic_cast<VANSOperator<dim, double> *>(
            this->system_operator.get()))
        {
          TimerOutput::Scope t(this->computing_timer,
                               "Prepare MF operator for VANS");

          mf_operator->compute_void_fraction(
            this->particle_projector.dof_handler,
            this->particle_projector.void_fraction_solution,
            this->time_derivative_void_fraction);

          mf_operator->compute_particle_fluid_interaction(
            this->particle_projector.fluid_force_on_particles_two_way_coupling
              .dof_handler,
            this->particle_projector.fluid_force_on_particles_two_way_coupling
              .particle_field_solution,
            this->particle_projector.fluid_drag_on_particles.dof_handler,
            this->particle_projector.fluid_drag_on_particles
              .particle_field_solution,
            this->particle_projector.particle_velocity.dof_handler,
            this->particle_projector.particle_velocity.particle_field_solution,
            this->particle_projector.momentum_transfer_coefficient.dof_handler,
            this->particle_projector.momentum_transfer_coefficient
              .particle_field_solution);
        }

      this->iterate();

      {
        announce_string(this->pcout, "DEM");
        TimerOutput::Scope t(this->computing_timer, "DEM_Iterator");

        // First DEM iteration of the CFD iteration
        dem_action_manager->first_dem_of_cfddem_iteration_step();

        // Load balancing if needed
        load_balance();

        // Update DEM time-step when the simulation uses adaptive time-stepping
        if (this->simulation_control->is_adaptive_time_stepping())
          update_dem_time_step();

        contact_search_counter = 0;
        for (unsigned int dem_counter = 0; dem_counter < coupling_frequency;
             ++dem_counter)
          {
            // dem_iterator carries out the particle-particle and
            // particle_wall force calculations, integration and
            // update_ghost
            dem_iterator(dem_counter);
          }
      }

      contact_search_total_number += contact_search_counter;

      this->pcout << "Finished " << coupling_frequency << " DEM iterations"
                  << std::endl;

      this->postprocess(false);
      this->postprocess_cfd_dem();
      this->finish_time_step_fd();

      if (this->cfd_dem_simulation_parameters.cfd_dem.particle_statistics)
        report_particle_statistics();
    }

  // Write particle-wall collision statistics file if enabled
  // if (dem_parameters.post_processing.particle_wall_collision_statistics)
  // write_collision_stats(dem_parameters,
  //                        collision_event_log,
  //                        this->mpi_communicator);
  if (this->simulation_parameters.timer.type == Parameters::Timer::Type::end)
    this->print_mg_setup_times();

  this->finish_simulation();
}

// Pre-compile the 2D and 3D CFD-DEM matrix-free solver
template class CFDDEMMatrixFree<2>;
template class CFDDEMMatrixFree<3>;
