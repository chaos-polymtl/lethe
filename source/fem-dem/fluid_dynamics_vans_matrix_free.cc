
// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>
#include <core/lethe_grid_tools.h>
#include <core/manifolds.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <dem/particle_handler_conversion.h>
#include <fem-dem/fluid_dynamics_vans_matrix_free.h>
#include <fem-dem/fluid_dynamics_vans_matrix_free_operators.h>

#include <deal.II/grid/grid_tools.h>

template <int dim>
MFNavierStokesVANSPreconditionGMG<dim>::MFNavierStokesVANSPreconditionGMG(
  const CFDDEMSimulationParameters<dim> &param,
  const DoFHandler<dim>                 &dof_handler,
  const DoFHandler<dim>                 &dof_handler_fe_q_iso_q1)
  : MFNavierStokesPreconditionGMG<dim>(param.cfd_parameters,
                                       dof_handler,
                                       dof_handler_fe_q_iso_q1)
  , cfd_dem_simulation_parameters(param)
{}

template <int dim>
void
MFNavierStokesVANSPreconditionGMG<dim>::create_level_operator(
  const unsigned int level)
{
  this->mg_operators[level] = std::make_shared<VANSOperator<dim, MGNumber>>(
    cfd_dem_simulation_parameters.cfd_dem);
}

template <int dim>
void
MFNavierStokesVANSPreconditionGMG<dim>::initialize(
  const std::shared_ptr<SimulationControl> &simulation_control,
  FlowControl<dim>                         &flow_control,
  const VectorType                         &present_solution,
  const VectorType                         &time_derivative_previous_solutions,
  const VectorType                         &time_derivative_void_fraction,
  const ParticleProjector<dim>             &particle_projector)
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::lsmg)
    {
      AssertThrow(false, ExcNotImplemented());
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::gcmg)
    {
      const unsigned int min_level = this->minlevel;
      const unsigned int max_level = this->maxlevel;

      // Resize all of the dof handlers related to every field for the
      // particle-fluid coupling
      this->void_fraction_dof_handlers.resize(min_level, max_level);
      this->pf_force_dof_handlers.resize(min_level, max_level);
      this->pf_drag_dof_handlers.resize(min_level, max_level);
      this->particle_velocity_dof_handlers.resize(min_level, max_level);
      this->momentum_transfer_coefficient_dof_handlers.resize(min_level,
                                                              max_level);


      for (unsigned int l = min_level; l <= max_level; l++)
        {
          // Void fraction
          this->void_fraction_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->void_fraction_dof_handlers[l].distribute_dofs(
            particle_projector.dof_handler.get_fe());

          // Particle-fluid force
          this->pf_force_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->pf_force_dof_handlers[l].distribute_dofs(
            particle_projector.fluid_force_on_particles_two_way_coupling
              .dof_handler.get_fe());

          // Particle-fluid drag
          this->pf_drag_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->pf_drag_dof_handlers[l].distribute_dofs(
            particle_projector.fluid_drag_on_particles.dof_handler.get_fe());

          // Particle velocity field
          this->particle_velocity_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->particle_velocity_dof_handlers[l].distribute_dofs(
            particle_projector.particle_velocity.dof_handler.get_fe());

          // Momentum exchange coefficient
          this->momentum_transfer_coefficient_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->momentum_transfer_coefficient_dof_handlers[l].distribute_dofs(
            particle_projector.momentum_transfer_coefficient.dof_handler
              .get_fe());
        }

      this->transfers_void_fraction.resize(min_level, max_level);
      this->transfers_pf_force.resize(min_level, max_level);
      this->transfers_pf_drag.resize(min_level, max_level);
      this->transfers_particle_velocity.resize(min_level, max_level);
      this->transfers_momentum_transfer_coefficient.resize(min_level,
                                                           max_level);


      for (unsigned int l = min_level; l < max_level; l++)
        {
          this->transfers_void_fraction[l + 1].reinit(
            this->void_fraction_dof_handlers[l + 1],
            this->void_fraction_dof_handlers[l],
            {},
            {});

          this->transfers_pf_force[l + 1].reinit(
            this->pf_force_dof_handlers[l + 1],
            this->pf_force_dof_handlers[l],
            {},
            {});

          this->transfers_pf_drag[l + 1].reinit(
            this->pf_drag_dof_handlers[l + 1],
            this->pf_drag_dof_handlers[l],
            {},
            {});

          this->transfers_particle_velocity[l + 1].reinit(
            this->particle_velocity_dof_handlers[l + 1],
            this->particle_velocity_dof_handlers[l],
            {},
            {});

          this->transfers_momentum_transfer_coefficient[l + 1].reinit(
            this->momentum_transfer_coefficient_dof_handlers[l + 1],
            this->momentum_transfer_coefficient_dof_handlers[l],
            {},
            {});
        }

      // Make the transfer operators for every field we will transfer
      this->mg_transfer_gc_void_fraction =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_void_fraction);

      this->mg_transfer_gc_pf_force =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_pf_force);

      this->mg_transfer_gc_pf_drag =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_pf_drag);

      this->mg_transfer_gc_particle_velocity =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_particle_velocity);

      this->mg_transfer_gc_momentum_transfer_coefficient =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_momentum_transfer_coefficient);


      // Build transfer operator for every field
      this->mg_transfer_gc_void_fraction->build(
        particle_projector.dof_handler, [&](const auto l, auto &vec) {
          vec.reinit(
            this->void_fraction_dof_handlers[l].locally_owned_dofs(),
            DoFTools::extract_locally_active_dofs(
              this->void_fraction_dof_handlers[l]),
            this->void_fraction_dof_handlers[l].get_mpi_communicator());
        });

      this->mg_transfer_gc_pf_force->build(
        particle_projector.fluid_force_on_particles_two_way_coupling
          .dof_handler,
        [&](const auto l, auto &vec) {
          vec.reinit(this->pf_force_dof_handlers[l].locally_owned_dofs(),
                     DoFTools::extract_locally_active_dofs(
                       this->pf_force_dof_handlers[l]),
                     this->pf_force_dof_handlers[l].get_mpi_communicator());
        });

      this->mg_transfer_gc_pf_force->build(
        particle_projector.fluid_drag_on_particles.dof_handler,
        [&](const auto l, auto &vec) {
          vec.reinit(this->pf_drag_dof_handlers[l].locally_owned_dofs(),
                     DoFTools::extract_locally_active_dofs(
                       this->pf_drag_dof_handlers[l]),
                     this->pf_drag_dof_handlers[l].get_mpi_communicator());
        });

      this->mg_transfer_gc_particle_velocity->build(
        particle_projector.particle_velocity.dof_handler,
        [&](const auto l, auto &vec) {
          vec.reinit(
            this->particle_velocity_dof_handlers[l].locally_owned_dofs(),
            DoFTools::extract_locally_active_dofs(
              this->particle_velocity_dof_handlers[l]),
            this->particle_velocity_dof_handlers[l].get_mpi_communicator());
        });

      this->mg_transfer_gc_momentum_transfer_coefficient->build(
        particle_projector.momentum_transfer_coefficient.dof_handler,
        [&](const auto l, auto &vec) {
          vec.reinit(this->momentum_transfer_coefficient_dof_handlers[l]
                       .locally_owned_dofs(),
                     DoFTools::extract_locally_active_dofs(
                       this->momentum_transfer_coefficient_dof_handlers[l]),
                     this->momentum_transfer_coefficient_dof_handlers[l]
                       .get_mpi_communicator());
        });

      // Create the MG Level Object for every field
      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_void_fraction_solution(this->minlevel, this->maxlevel);

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_time_derivative_void_fraction(this->minlevel, this->maxlevel);

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_pf_forces_solution(this->minlevel, this->maxlevel);

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_pf_drag_solution(this->minlevel, this->maxlevel);

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_particle_velocity_solution(this->minlevel, this->maxlevel);

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_momentum_transfer_coefficient_solution(this->minlevel,
                                                  this->maxlevel);

      // A deal.II vector is required here, so we take the deal.II vector
      // solution from the particle projector instead of the Trillinos vector
      // one. This is done for every field.

      // Void fraction
      this->mg_transfer_gc_void_fraction->interpolate_to_mg(
        particle_projector.dof_handler,
        mg_void_fraction_solution,
        particle_projector.void_fraction_solution);

      this->mg_transfer_gc_void_fraction->interpolate_to_mg(
        particle_projector.dof_handler,
        mg_time_derivative_void_fraction,
        time_derivative_void_fraction);

      particle_projector.fluid_force_on_particles_two_way_coupling
        .particle_field_solution.update_ghost_values();

      // Particle-fluid force
      this->mg_transfer_gc_pf_force->interpolate_to_mg(
        particle_projector.fluid_force_on_particles_two_way_coupling
          .dof_handler,
        mg_pf_forces_solution,
        particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution);

      // Particle-fluid drag
      particle_projector.fluid_drag_on_particles.particle_field_solution
        .update_ghost_values();

      this->mg_transfer_gc_pf_force->interpolate_to_mg(
        particle_projector.fluid_drag_on_particles.dof_handler,
        mg_pf_drag_solution,
        particle_projector.fluid_drag_on_particles.particle_field_solution);

      // Particle-velocity
      particle_projector.particle_velocity.particle_field_solution
        .update_ghost_values();

      this->mg_transfer_gc_pf_force->interpolate_to_mg(
        particle_projector.particle_velocity.dof_handler,
        mg_particle_velocity_solution,
        particle_projector.particle_velocity.particle_field_solution);

      // Momentum transfer coefficient
      particle_projector.momentum_transfer_coefficient.particle_field_solution
        .update_ghost_values();

      this->mg_transfer_gc_momentum_transfer_coefficient->interpolate_to_mg(
        particle_projector.momentum_transfer_coefficient.dof_handler,
        mg_momentum_transfer_coefficient_solution,
        particle_projector.momentum_transfer_coefficient
          .particle_field_solution);


      for (unsigned int l = min_level; l <= max_level; l++)
        {
          mg_void_fraction_solution[l].update_ghost_values();
          mg_time_derivative_void_fraction[l].update_ghost_values();
          mg_pf_forces_solution[l].update_ghost_values();
          mg_pf_drag_solution[l].update_ghost_values();
          mg_particle_velocity_solution[l].update_ghost_values();
          mg_momentum_transfer_coefficient_solution[l].update_ghost_values();

          if (auto mf_operator = dynamic_cast<VANSOperator<dim, double> *>(
                &(*this->mg_operators[l])))
            {
              mf_operator->compute_void_fraction(
                this->void_fraction_dof_handlers[l],
                mg_void_fraction_solution[l],
                mg_time_derivative_void_fraction[l]);

              mf_operator->compute_particle_fluid_interaction(
                this->pf_force_dof_handlers[l],
                mg_pf_forces_solution[l],
                this->pf_drag_dof_handlers[l],
                mg_pf_drag_solution[l],
                this->particle_velocity_dof_handlers[l],
                mg_particle_velocity_solution[l],
                this->momentum_transfer_coefficient_dof_handlers[l],
                mg_momentum_transfer_coefficient_solution[l]);
            }
        }
    }

  MFNavierStokesPreconditionGMG<dim>::initialize(
    simulation_control,
    flow_control,
    present_solution,
    time_derivative_previous_solutions);
}

template <int dim>
FluidDynamicsVANSMatrixFree<dim>::FluidDynamicsVANSMatrixFree(
  CFDDEMSimulationParameters<dim> &param)
  : FluidDynamicsMatrixFree<dim>(param.cfd_parameters)
  , cfd_dem_simulation_parameters(param)
  , particle_mapping(1)
  , particle_handler(*this->triangulation,
                     particle_mapping,
                     DEM::CFDDEMProperties::n_properties)
  , particle_projector(
      &(*this->triangulation),
      param.void_fraction,
      this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver.at(
        PhysicsID::void_fraction),
      &particle_handler,
      this->cfd_dem_simulation_parameters.cfd_parameters.fem_parameters
        .void_fraction_order,
      this->cfd_dem_simulation_parameters.cfd_parameters.mesh.simplex,
      this->pcout)
  , has_periodic_boundaries(false)
{
  unsigned int n_pbc = 0;
  for (auto const &[id, type] :
       cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::periodic)
        {
          if (n_pbc++ > 1)
            {
              throw std::runtime_error(
                "GLS VANS solver does not support more than one periodic boundary condition.");
            }
          else
            {
              has_periodic_boundaries = true;
            }
        }
    }

  // The default MatrixFree solver sets a system_operator. We override the
  // Navier-Stokes operator with the volume-averaged Navier-Stokes operator.
  this->system_operator = std::make_shared<VANSOperator<dim, double>>(
    cfd_dem_simulation_parameters.cfd_dem);
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::setup_dofs()
{
  FluidDynamicsMatrixFree<dim>::setup_dofs();

  particle_projector.setup_dofs();
  particle_projector.setup_constraints(
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  // Initialize the time derivative of the void fraction.
  this->time_derivative_void_fraction.reinit(
    particle_projector.dof_handler.locally_owned_dofs(),
    DoFTools::extract_locally_relevant_dofs(particle_projector.dof_handler),
    this->mpi_communicator);
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::read_dem()
{
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
    *this->triangulation, particle_mapping, DEM::DEMProperties::n_properties);

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
          this->triangulation->load(filename.c_str());

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
      // function which requires a pararallel::distributed::triangulation
      convert_particle_handler<dim,
                               DEM::DEMProperties::PropertiesIndex,
                               DEM::CFDDEMProperties::PropertiesIndex>(
        *parallel_triangulation, temporary_particle_handler, particle_handler);

      // Exchange the ghost particles so that the particle handler is ready to
      // be used to calculate the initial conditions of the void fraction or
      // anything that requires the ghost particles.
      particle_handler.exchange_ghost_particles(true);
    }
  else
    {
      throw std::runtime_error(
        "The VANS application currently does not support "
        "triangulations other than parallel::distributed");
    }
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::assemble_system_rhs()
{
  // First we update the particle-fluid coupling terms.
  // The base matrix-free operator is not aware of the various VANS
  // coupling terms. We must do a cast here to ensure that the operator is
  // of the right type.
  if (auto mf_operator =
        dynamic_cast<VANSOperator<dim, double> *>(this->system_operator.get()))
    {
      TimerOutput::Scope t(this->computing_timer,
                           "Prepare MF operator for VANS");

      // If the coupling is implicit, we need to carry out the projection of
      // the particle information onto the CFD mesh at every RHS evaluation
      // since the momentum exchange term evolves.
      if (cfd_dem_simulation_parameters.cfd_dem.drag_coupling ==
          Parameters::DragCoupling::fully_implicit)
        particle_projector.calculate_particle_fluid_forces_projection(
          this->cfd_dem_simulation_parameters.cfd_dem,
          *this->dof_handler,
          this->evaluation_point,
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

      mf_operator->compute_particle_fluid_interaction(
        particle_projector.fluid_force_on_particles_two_way_coupling
          .dof_handler,
        particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution,
        particle_projector.fluid_drag_on_particles.dof_handler,
        particle_projector.fluid_drag_on_particles.particle_field_solution,
        particle_projector.particle_velocity.dof_handler,
        particle_projector.particle_velocity.particle_field_solution,
        particle_projector.momentum_transfer_coefficient.dof_handler,
        particle_projector.momentum_transfer_coefficient
          .particle_field_solution);
    }

  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");


  // Update the precomputed values needed for the evaluation of the residual.
  // This is needed, otherwise the line-search mechanism used in the Newton
  // method might fail even though the Newton step should have been accepted
  // due to a wrong evaluation of the residual and, consequently, a wrong
  // evaluation of the step length.
  this->evaluation_point.update_ghost_values();
  this->system_operator->evaluate_non_linear_term_and_calculate_tau(
    this->evaluation_point);

  this->system_operator->evaluate_residual(this->system_rhs,
                                           this->evaluation_point);

  // If mortar is enabled, add RHS entries
  if (this->simulation_parameters.mortar_parameters.enable)
    {
      this->system_operator->mortar_coupling_operator_mf->vmult_add(
        this->system_rhs, this->evaluation_point);
      this->system_rhs.compress(VectorOperation::add);
    }

  this->system_rhs *= -1.0;

  // Provide residual to simulation control for stopping criterion when using
  // steady bdf
  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::finish_time_step_fd()
{
  // Void fraction percolation must be done before the time step is finished to
  // ensure that the checkpointed information is correct
  particle_projector.percolate_void_fraction();

  FluidDynamicsMatrixFree<dim>::finish_time_step();
}

template <int dim>
std::vector<OutputStruct<dim, LinearAlgebra::distributed::Vector<double>>>
FluidDynamicsVANSMatrixFree<dim>::gather_output_hook()
{
  // Make sure all of the particle-related vectors have adequate ghost values.
  particle_projector.void_fraction_solution.update_ghost_values();
  particle_projector.fluid_drag_on_particles.particle_field_solution
    .update_ghost_values();
  particle_projector.particle_velocity.particle_field_solution
    .update_ghost_values();
  particle_projector.momentum_transfer_coefficient.particle_field_solution
    .update_ghost_values();

  std::vector<std::string> name = {"void_fraction"};
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation = {
      DataComponentInterpretation::component_is_scalar};
  OutputStructSolution<dim, LinearAlgebra::distributed::Vector<double>>
    void_fraction_struct(particle_projector.dof_handler,
                         particle_projector.void_fraction_solution,
                         name,
                         component_interpretation);

  std::vector<std::string> force_names(dim, "fluid_drag_on_particles");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    force_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  OutputStructSolution<dim, LinearAlgebra::distributed::Vector<double>>
    particle_fluid_drag_struct(
      particle_projector.fluid_drag_on_particles.dof_handler,
      particle_projector.fluid_drag_on_particles.particle_field_solution,
      force_names,
      force_component_interpretation);

  std::vector<std::string> particle_velocity_names(dim, "particle_velocity");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    particle_velocity_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);

  OutputStructSolution<dim, LinearAlgebra::distributed::Vector<double>>
    particle_velocity_struct(
      particle_projector.particle_velocity.dof_handler,
      particle_projector.particle_velocity.particle_field_solution,
      particle_velocity_names,
      particle_velocity_component_interpretation);

  std::vector<std::string> mtc_name(1, "momentum_transfer_coefficient");
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    mtc_component_interpretation(
      1, DataComponentInterpretation::component_is_scalar);

  OutputStructSolution<dim, LinearAlgebra::distributed::Vector<double>>
    mtc_struct(
      particle_projector.momentum_transfer_coefficient.dof_handler,
      particle_projector.momentum_transfer_coefficient.particle_field_solution,
      mtc_name,
      mtc_component_interpretation);


  return {void_fraction_struct,
          particle_fluid_drag_struct,
          particle_velocity_struct,
          mtc_struct};
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::evaluate_time_derivative_void_fraction()
{
  this->time_derivative_void_fraction = 0;

  // If the void fraction time derivative is disabled, we just leave it at zero.
  if (cfd_dem_simulation_parameters.cfd_dem.void_fraction_time_derivative ==
      false)
    {
      // Even though it is zero, we still need to update its ghost values.
      time_derivative_void_fraction.update_ghost_values();
      return;
    }
  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();

  this->time_derivative_void_fraction.add(
    bdf_coefs[0], this->particle_projector.void_fraction_solution);
  for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
    {
      this->time_derivative_void_fraction.add(
        bdf_coefs[p + 1],
        this->particle_projector.void_fraction_previous_solution[p]);
    }

  // After the void fraction time derivative is calculated, we update its ghost
  // values.
  time_derivative_void_fraction.update_ghost_values();
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::solve()
{
  this->computing_timer.enter_subsection("Read mesh, manifolds and particles");

  read_mesh_and_manifolds(
    *this->triangulation,
    this->cfd_dem_simulation_parameters.cfd_parameters.mesh,
    this->cfd_dem_simulation_parameters.cfd_parameters.manifolds_parameters,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
        .restart ||
      this->cfd_dem_simulation_parameters.void_fraction->read_dem == true,
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);

  if (this->cfd_dem_simulation_parameters.void_fraction->read_dem == true &&
      this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
          .restart == false)
    read_dem();

  this->computing_timer.leave_subsection("Read mesh, manifolds and particles");

  this->setup_dofs();

  particle_projector.initialize_void_fraction(
    this->simulation_control->get_current_time());

  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);

  // Only needed if other physics apart from fluid dynamics are enabled.
  if (this->multiphysics->get_active_physics().size() > 1)
    this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      if (this->forcing_function)
        this->forcing_function->set_time(
          this->simulation_control->get_current_time());

      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (!this->simulation_control->is_at_start())
        {
          this->refine_mesh();
        }

      if (time_stepping_is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->computing_timer.enter_subsection(
            "Calculate time derivative previous solutions");

          this->calculate_time_derivative_previous_solutions();
          this->time_derivative_previous_solutions.update_ghost_values();
          this->system_operator->evaluate_time_derivative_previous_solutions(
            this->time_derivative_previous_solutions);

          this->computing_timer.leave_subsection(
            "Calculate time derivative previous solutions");

          if (this->simulation_parameters.flow_control.enable_flow_control)
            this->system_operator->update_beta_force(
              this->flow_control.get_beta());
        }

      {
        TimerOutput::Scope t(this->computing_timer,
                             "Calculate particle-fluid projection");

        // Calculate the new void fraction
        particle_projector.calculate_void_fraction(
          this->simulation_control->get_current_time());

        // and its derivative
        evaluate_time_derivative_void_fraction();

        particle_projector.calculate_particle_fluid_forces_projection(
          this->cfd_dem_simulation_parameters.cfd_dem,
          *this->dof_handler,
          *this->present_solution,
          *this->previous_solutions,
          cfd_dem_simulation_parameters.dem_parameters
            .lagrangian_physical_properties.g,
          NavierStokesScratchData<dim>(
            this->simulation_control,
            this->simulation_parameters.physical_properties_manager,
            *this->fe,
            *this->cell_quadrature,
            *this->mapping,
            *this->face_quadrature));

        // The base matrix-free operator is not aware of the various VANS
        // coupling term. We must do a cast here to ensure that the operator is
        // of the right type.
        if (auto mf_operator = dynamic_cast<VANSOperator<dim, double> *>(
              this->system_operator.get()))
          {
            TimerOutput::Scope t(this->computing_timer,
                                 "Prepare MF operator for VANS");

            mf_operator->compute_void_fraction(
              particle_projector.dof_handler,
              particle_projector.void_fraction_solution,
              this->time_derivative_void_fraction);

            mf_operator->compute_particle_fluid_interaction(
              particle_projector.fluid_force_on_particles_two_way_coupling
                .dof_handler,
              particle_projector.fluid_force_on_particles_two_way_coupling
                .particle_field_solution,
              particle_projector.fluid_drag_on_particles.dof_handler,
              particle_projector.fluid_drag_on_particles
                .particle_field_solution,
              particle_projector.particle_velocity.dof_handler,
              particle_projector.particle_velocity.particle_field_solution,
              particle_projector.momentum_transfer_coefficient.dof_handler,
              particle_projector.momentum_transfer_coefficient
                .particle_field_solution);
          }
      }

      this->iterate();
      this->postprocess(false);
      this->finish_time_step();

      if (this->simulation_parameters.timer.type ==
          Parameters::Timer::Type::iteration)
        this->print_mg_setup_times();
    }
  if (this->simulation_parameters.timer.type == Parameters::Timer::Type::end)
    this->print_mg_setup_times();

  this->finish_simulation();
}


template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::create_GMG()
{
  this->gmg_preconditioner =
    std::make_shared<MFNavierStokesVANSPreconditionGMG<dim>>(
      this->cfd_dem_simulation_parameters,
      *this->dof_handler,
      this->dof_handler_fe_q_iso_q1);

  this->gmg_preconditioner->reinit(this->mapping,
                                   this->cell_quadrature,
                                   this->forcing_function,
                                   this->simulation_control,
                                   this->physical_properties_manager,
                                   this->fe);
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::initialize_GMG()
{
  dynamic_cast<MFNavierStokesVANSPreconditionGMG<dim> *>(
    this->gmg_preconditioner.get())
    ->initialize(this->simulation_control,
                 this->flow_control,
                 *this->present_solution,
                 this->time_derivative_previous_solutions,
                 this->time_derivative_void_fraction,
                 this->particle_projector);
}

// Pre-compile the 2D and 3D solver to ensure that the
// library is valid before we actually compile the solver
template class FluidDynamicsVANSMatrixFree<2>;
template class FluidDynamicsVANSMatrixFree<3>;
