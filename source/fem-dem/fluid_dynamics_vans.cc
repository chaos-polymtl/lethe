// SPDX-FileCopyrightText: Copyright (c) 2020-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>
#include <core/lethe_grid_tools.h>

#include <dem/particle_handler_conversion.h>
#include <fem-dem/fluid_dynamics_vans.h>
#include <fem-dem/particle_projector.h>

#include <deal.II/base/work_stream.h>

// Constructor for class FluidDynamicsVANS
template <int dim>
FluidDynamicsVANS<dim>::FluidDynamicsVANS(
  CFDDEMSimulationParameters<dim> &nsparam)
  : FluidDynamicsMatrixBased<dim>(nsparam.cfd_parameters)
  , cfd_dem_simulation_parameters(nsparam)
  , particle_mapping(1)
  , particle_handler(*this->triangulation,
                     particle_mapping,
                     DEM::CFDDEMProperties::n_properties)
  , particle_projector(
      &(*this->triangulation),
      nsparam.void_fraction,
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
}

template <int dim>
FluidDynamicsVANS<dim>::~FluidDynamicsVANS()
{
  this->dof_handler->clear();
}

template <int dim>
void
FluidDynamicsVANS<dim>::setup_dofs()
{
  FluidDynamicsMatrixBased<dim>::setup_dofs();

  particle_projector.setup_dofs();
  particle_projector.setup_constraints(
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);
}

template <int dim>
void
FluidDynamicsVANS<dim>::finish_time_step_fd()
{
  // Void fraction percolation must be done before the time step is finished to
  // ensure that the checkpointed information is correct
  particle_projector.percolate_void_fraction();

  FluidDynamicsMatrixBased<dim>::finish_time_step();
}

template <int dim>
void
FluidDynamicsVANS<dim>::read_dem()
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
    }
  else
    {
      throw std::runtime_error(
        "VANS equations currently do not support "
        "triangulations other than parallel::distributed");
    }
}

template <int dim>
void
FluidDynamicsVANS<dim>::calculate_void_fraction(const double time)
{
  TimerOutput::Scope t(this->computing_timer, "Calculate void fraction");
  particle_projector.calculate_void_fraction(time);
}

template <int dim>
void
FluidDynamicsVANS<dim>::vertices_cell_mapping()
{
  // Find all the cells around each vertex
  TimerOutput::Scope t(this->computing_timer, "Map vertices to cell");

  LetheGridTools::vertices_cell_mapping(this->particle_projector.dof_handler,
                                        vertices_to_cell);

  if (has_periodic_boundaries)
    LetheGridTools::vertices_cell_mapping_with_periodic_boundaries(
      this->particle_projector.dof_handler, vertices_to_periodic_cell);
}

// Do an iteration with the NavierStokes Solver
// Handles the fact that we may or may not be at a first
// iteration with the solver and sets the initial conditions
template <int dim>
void
FluidDynamicsVANS<dim>::iterate()
{
  announce_string(this->pcout, "Volume-Averaged Fluid Dynamics");
  this->forcing_function->set_time(
    this->simulation_control->get_current_time());

  PhysicsSolver<GlobalVectorType>::solve_governing_system();
}

template <int dim>
void
FluidDynamicsVANS<dim>::setup_assemblers()
{
  this->assemblers.clear();
  particle_fluid_assemblers.clear();

  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::function_weak))
    {
      this->assemblers.push_back(
        std::make_shared<WeakDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(
        BoundaryConditions::BoundaryType::partial_slip))
    {
      this->assemblers.push_back(
        std::make_shared<PartialSlipDirichletBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::outlet))
    {
      this->assemblers.push_back(std::make_shared<OutletBoundaryCondition<dim>>(
        this->simulation_control,
        this->simulation_parameters.boundary_conditions));
    }
  if (this->check_existance_of_bc(BoundaryConditions::BoundaryType::pressure))
    {
      this->assemblers.push_back(
        std::make_shared<PressureBoundaryCondition<dim>>(
          this->simulation_control,
          this->simulation_parameters.boundary_conditions));
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.drag_force == true)
    {
      // Particle_Fluid Interactions Assembler
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::difelice)
        {
          // DiFelice Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerDiFelice<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::rong)
        {
          // Rong Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerRong<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::dallavalle)
        {
          // Dallavalle Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerDallavalle<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }

      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::kochhill)
        {
          // Koch and Hill Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerKochHill<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::beetstra)
        {
          // Beetstra drag model assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerBeetstra<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_model ==
          Parameters::DragModel::gidaspow)
        {
          // Gidaspow Model drag Assembler
          particle_fluid_assemblers.push_back(
            std::make_shared<VANSAssemblerGidaspow<dim>>(
              this->cfd_dem_simulation_parameters.cfd_dem));
        }
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.saffman_lift_force == true)
    // Saffman Mei Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerSaffmanMei<dim>>());

  if (this->cfd_dem_simulation_parameters.cfd_dem.magnus_lift_force == true)
    // Magnus Lift Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerMagnus<dim>>());

  if (this->cfd_dem_simulation_parameters.cfd_dem.rotational_viscous_torque ==
      true)
    // Viscous Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerViscousTorque<dim>>());

  if (this->cfd_dem_simulation_parameters.cfd_dem.vortical_viscous_torque ==
      true)
    // Vortical Torque Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerVorticalTorque<dim>>());

  if (this->cfd_dem_simulation_parameters.cfd_dem.buoyancy_force == true)
    // Buoyancy Force Assembler
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerBuoyancy<dim>>(
        this->cfd_dem_simulation_parameters.dem_parameters
          .lagrangian_physical_properties.g));

  if (this->cfd_dem_simulation_parameters.cfd_dem.pressure_force == true)
    // Pressure Force
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerPressureForce<dim>>(
        this->cfd_dem_simulation_parameters.cfd_dem));

  if (this->cfd_dem_simulation_parameters.cfd_dem.shear_force == true)
    // Shear Force
    particle_fluid_assemblers.push_back(
      std::make_shared<VANSAssemblerShearForce<dim>>(
        this->cfd_dem_simulation_parameters.cfd_dem));

  // Time-stepping schemes
  if (time_stepping_is_bdf(this->simulation_control->get_assembly_method()))
    {
      this->assemblers.push_back(std::make_shared<VANSAssemblerBDF<dim>>(
        this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));
    }

  //  Fluid_Particle Interactions Assembler
  if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
    {
      this->assemblers.push_back(
        std::make_shared<VANSAssemblerFPIProjection<dim>>(
          this->cfd_dem_simulation_parameters.cfd_dem));
    }
  else
    {
      this->assemblers.push_back(std::make_shared<VANSAssemblerFPI<dim>>(
        this->cfd_dem_simulation_parameters.cfd_dem));
    }

  // The core assembler should always be the last assembler to be called
  // in the stabilized formulation as to have all strong residual and
  // jacobian stored. Core assembler
  if (this->cfd_dem_simulation_parameters.cfd_dem.vans_model ==
      Parameters::VANSModel::modelA)
    this->assemblers.push_back(std::make_shared<VANSAssemblerCoreModelA<dim>>(
      this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));

  if (this->cfd_dem_simulation_parameters.cfd_dem.vans_model ==
      Parameters::VANSModel::modelB)
    this->assemblers.push_back(std::make_shared<VANSAssemblerCoreModelB<dim>>(
      this->simulation_control, this->cfd_dem_simulation_parameters.cfd_dem));
}

template <int dim>
void
FluidDynamicsVANS<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
  this->system_matrix = 0;

  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_control,
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> &dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof.get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          const DoFHandler<dim> &projected_phase_fraction_gradient_dof_handler =
            this->multiphysics
              ->get_projected_phase_fraction_gradient_dof_handler();
          const DoFHandler<dim> &curvature_dof_handler =
            this->multiphysics->get_curvature_dof_handler();
          scratch_data.enable_projected_phase_fraction_gradient(
            projected_phase_fraction_gradient_dof_handler.get_fe(),
            *this->cell_quadrature,
            *this->mapping);
          scratch_data.enable_curvature(curvature_dof_handler.get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
        }
    }

  scratch_data.enable_void_fraction(*particle_projector.fe,
                                    *this->cell_quadrature,
                                    *this->mapping);
  if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
    {
      scratch_data.enable_particle_field_projection(
        *this->cell_quadrature,
        *this->mapping,
        *particle_projector.fluid_drag_on_particles.fe,
        *particle_projector.fluid_force_on_particles_two_way_coupling.fe,
        *particle_projector.particle_velocity.fe,
        *particle_projector.momentum_transfer_coefficient.fe);
    }
  scratch_data.enable_particle_fluid_interactions(
    particle_handler.n_global_max_particles_per_cell(),
    this->cfd_dem_simulation_parameters.cfd_dem.interpolated_void_fraction);

  WorkStream::run(
    this->dof_handler->begin_active(),
    this->dof_handler->end(),
    *this,
    &FluidDynamicsVANS::assemble_local_system_matrix,
    &FluidDynamicsVANS::copy_local_matrix_to_global_matrix,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));
  this->system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
FluidDynamicsVANS<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    *this->previous_solutions,
    this->sdirk_vectors.sum_over_previous_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  std::optional<typename DoFHandler<dim>::active_cell_iterator> phase_cell;

  if (this->simulation_parameters.multiphysics.VOF)
    {
      phase_cell.emplace(&(*(this->triangulation)),
                         cell->level(),
                         cell->index(),
                         &this->multiphysics->get_dof_handler(PhysicsID::VOF));

      scratch_data.reinit_vof(
        *phase_cell,
        this->multiphysics->get_solution(PhysicsID::VOF),
        this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->particle_projector.dof_handler);

  scratch_data.reinit_void_fraction(
    void_fraction_cell,
    particle_projector.void_fraction_locally_relevant,
    particle_projector.previous_void_fraction);

  scratch_data.calculate_physical_properties();

  if (this->simulation_parameters.multiphysics.VOF)
    {
      scratch_data.reinit_particle_fluid_interactions(
        cell,
        void_fraction_cell,
        *phase_cell,
        this->evaluation_point,
        (*this->previous_solutions)[0],
        this->particle_projector.void_fraction_locally_relevant,
        particle_handler,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling,
        this->multiphysics->get_filtered_solution(PhysicsID::VOF));
    }
  else
    {
      scratch_data.reinit_particle_fluid_interactions(
        cell,
        void_fraction_cell,
        this->evaluation_point,
        (*this->previous_solutions)[0],
        this->particle_projector.void_fraction_locally_relevant,
        particle_handler,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling);
    }

  for (auto &pf_assembler : particle_fluid_assemblers)
    {
      pf_assembler->calculate_particle_fluid_interactions(scratch_data);
    }
  if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
    {
      typename DoFHandler<dim>::active_cell_iterator particle_drag_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &this->particle_projector.fluid_drag_on_particles.dof_handler);

      typename DoFHandler<dim>::active_cell_iterator
        particle_two_way_coupling_force_cell(
          &(*(this->triangulation)),
          cell->level(),
          cell->index(),
          &this->particle_projector.fluid_force_on_particles_two_way_coupling
             .dof_handler);

      typename DoFHandler<dim>::active_cell_iterator particle_velocity_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &this->particle_projector.particle_velocity.dof_handler);

      typename DoFHandler<dim>::active_cell_iterator
        particle_momentum_transfer_coefficient_cell(
          &(*(this->triangulation)),
          cell->level(),
          cell->index(),
          &this->particle_projector.momentum_transfer_coefficient.dof_handler);

      scratch_data.calculate_particle_fields_values(
        particle_drag_cell,
        particle_two_way_coupling_force_cell,
        particle_velocity_cell,
        particle_momentum_transfer_coefficient_cell,
        particle_projector.fluid_drag_on_particles.particle_field_solution,
        particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution,
        particle_projector.particle_velocity.particle_field_solution,
        particle_projector.momentum_transfer_coefficient
          .particle_field_solution,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling);
    }

  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
FluidDynamicsVANS<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              this->system_matrix);
}

template <int dim>
void
FluidDynamicsVANS<dim>::assemble_system_rhs()
{
  this->computing_timer.enter_subsection("Assemble RHS");
  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_control,
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);
  this->computing_timer.leave_subsection("Assemble RHS");

  if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
    {
      scratch_data.enable_particle_field_projection(
        *this->cell_quadrature,
        *this->mapping,
        *particle_projector.fluid_drag_on_particles.fe,
        *particle_projector.fluid_force_on_particles_two_way_coupling.fe,
        *particle_projector.particle_velocity.fe,
        *particle_projector.momentum_transfer_coefficient.fe);
      if (this->cfd_dem_simulation_parameters.cfd_dem.drag_coupling ==
          Parameters::DragCoupling::fully_implicit)
        {
          TimerOutput::Scope t(this->computing_timer,
                               "Calculate particle-fluid projection");
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
        }
    }

  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  this->system_rhs = 0;

  setup_assemblers();

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> &dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof.get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);
    }
  scratch_data.enable_void_fraction(*particle_projector.fe,
                                    *this->cell_quadrature,
                                    *this->mapping);
  scratch_data.enable_particle_fluid_interactions(
    particle_handler.n_global_max_particles_per_cell(),
    this->cfd_dem_simulation_parameters.cfd_dem.interpolated_void_fraction);

  WorkStream::run(
    this->dof_handler->begin_active(),
    this->dof_handler->end(),
    *this,
    &FluidDynamicsVANS::assemble_local_system_rhs,
    &FluidDynamicsVANS::copy_local_rhs_to_global_rhs,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
FluidDynamicsVANS<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim>                         &scratch_data,
  StabilizedMethodsTensorCopyData<dim>                 &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    *this->previous_solutions,
    this->sdirk_vectors.sum_over_previous_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  std::optional<typename DoFHandler<dim>::active_cell_iterator> phase_cell;

  if (this->simulation_parameters.multiphysics.VOF)
    {
      phase_cell.emplace(&(*(this->triangulation)),
                         cell->level(),
                         cell->index(),
                         &this->multiphysics->get_dof_handler(PhysicsID::VOF));

      scratch_data.reinit_vof(
        *phase_cell,
        this->multiphysics->get_solution(PhysicsID::VOF),
        this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
    &(*(this->triangulation)),
    cell->level(),
    cell->index(),
    &this->particle_projector.dof_handler);

  scratch_data.reinit_void_fraction(
    void_fraction_cell,
    particle_projector.void_fraction_locally_relevant,
    particle_projector.previous_void_fraction);

  scratch_data.calculate_physical_properties();

  if (this->simulation_parameters.multiphysics.VOF)
    {
      scratch_data.reinit_particle_fluid_interactions(
        cell,
        void_fraction_cell,
        *phase_cell,
        this->evaluation_point,
        (*this->previous_solutions)[0],
        particle_projector.void_fraction_locally_relevant,
        particle_handler,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling,
        this->multiphysics->get_filtered_solution(PhysicsID::VOF));
    }
  else
    {
      scratch_data.reinit_particle_fluid_interactions(
        cell,
        void_fraction_cell,
        this->evaluation_point,
        (*this->previous_solutions)[0],
        particle_projector.void_fraction_locally_relevant,
        particle_handler,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling);
    }
  for (auto &pf_assembler : particle_fluid_assemblers)
    {
      pf_assembler->calculate_particle_fluid_interactions(scratch_data);
    }

  if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
    {
      typename DoFHandler<dim>::active_cell_iterator particle_drag_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &this->particle_projector.fluid_drag_on_particles.dof_handler);

      typename DoFHandler<dim>::active_cell_iterator
        particle_two_way_coupling_force_cell(
          &(*(this->triangulation)),
          cell->level(),
          cell->index(),
          &this->particle_projector.fluid_force_on_particles_two_way_coupling
             .dof_handler);

      typename DoFHandler<dim>::active_cell_iterator particle_velocity_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &this->particle_projector.particle_velocity.dof_handler);

      typename DoFHandler<dim>::active_cell_iterator
        particle_momentum_transfer_coefficient_cell(
          &(*(this->triangulation)),
          cell->level(),
          cell->index(),
          &this->particle_projector.momentum_transfer_coefficient.dof_handler);

      scratch_data.calculate_particle_fields_values(
        particle_drag_cell,
        particle_two_way_coupling_force_cell,
        particle_velocity_cell,
        particle_momentum_transfer_coefficient_cell,
        particle_projector.fluid_drag_on_particles.particle_field_solution,
        particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution,
        particle_projector.particle_velocity.particle_field_solution,
        particle_projector.momentum_transfer_coefficient
          .particle_field_solution,
        cfd_dem_simulation_parameters.cfd_dem.drag_coupling);
    }

  copy_data.reset();

  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
FluidDynamicsVANS<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

template <int dim>
std::vector<OutputStruct<dim, GlobalVectorType>>
FluidDynamicsVANS<dim>::gather_output_hook()
{
  std::vector<OutputStruct<dim, GlobalVectorType>> solution_output_structs;
  solution_output_structs.emplace_back(
    std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
    particle_projector.dof_handler,
    particle_projector.void_fraction_locally_relevant,
    std::vector<std::string>{"void_fraction"},
    std::vector<DataComponentInterpretation::DataComponentInterpretation>{
      DataComponentInterpretation::component_is_scalar});
  if (this->cfd_dem_simulation_parameters.void_fraction
        ->project_particle_velocity)
    {
      std::vector<std::string> names(dim, "particle_velocity");
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);

#ifndef LETHE_USE_LDV
      // Since the particle velocity field is now only a deal.II distributed
      // vector, we create a temporary GlobalVectorType (a Trilinos vector) and
      // copy the content into it.
      GlobalVectorType particle_velocity;
      particle_velocity.reinit(
        particle_projector.particle_velocity.locally_owned_dofs,
        particle_projector.particle_velocity.locally_relevant_dofs,
        this->mpi_communicator);

      convert_vector_dealii_to_trilinos(
        particle_velocity,
        particle_projector.particle_velocity.particle_field_solution);
      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        particle_projector.particle_velocity.dof_handler,
        particle_velocity,
        names,
        data_interpretation);

      if (this->cfd_dem_simulation_parameters.cfd_dem.project_particle_forces)
        {
          GlobalVectorType fluid_drag_on_particles;
          fluid_drag_on_particles.reinit(
            particle_projector.particle_velocity.locally_owned_dofs,
            particle_projector.particle_velocity.locally_relevant_dofs,
            this->mpi_communicator);

          convert_vector_dealii_to_trilinos(
            fluid_drag_on_particles,
            particle_projector.fluid_drag_on_particles.particle_field_solution);

          solution_output_structs.emplace_back(
            std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
            this->particle_projector.fluid_drag_on_particles.dof_handler,
            fluid_drag_on_particles,
            std::vector<std::string>(dim, "Particle_drag"),
            std::vector<
              DataComponentInterpretation::DataComponentInterpretation>(
              dim, DataComponentInterpretation::component_is_part_of_vector));

          GlobalVectorType fluid_force_on_particles_two_way_coupling;
          fluid_force_on_particles_two_way_coupling.reinit(
            particle_projector.particle_velocity.locally_owned_dofs,
            particle_projector.particle_velocity.locally_relevant_dofs,
            this->mpi_communicator);

          convert_vector_dealii_to_trilinos(
            fluid_force_on_particles_two_way_coupling,
            particle_projector.fluid_force_on_particles_two_way_coupling
              .particle_field_solution);

          solution_output_structs.emplace_back(
            std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
            this->particle_projector.fluid_force_on_particles_two_way_coupling
              .dof_handler,
            fluid_force_on_particles_two_way_coupling,
            std::vector<std::string>(dim, "Particle_two_way_coupling_force"),
            std::vector<
              DataComponentInterpretation::DataComponentInterpretation>(
              dim, DataComponentInterpretation::component_is_part_of_vector));
        }
#else
      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        particle_projector.particle_velocity.dof_handler,
        particle_projector.particle_velocity.particle_field_solution,
        names,
        data_interpretation);

      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        this->particle_projector.fluid_drag_on_particles.dof_handler,
        this->particle_projector.fluid_drag_on_particles
          .particle_field_solution,
        std::vector<std::string>(dim, "Particle_drag"),
        std::vector<DataComponentInterpretation::DataComponentInterpretation>(
          dim, DataComponentInterpretation::component_is_part_of_vector));

      solution_output_structs.emplace_back(
        std::in_place_type<OutputStructSolution<dim, GlobalVectorType>>,
        this->particle_projector.fluid_force_on_particles_two_way_coupling
          .dof_handler,
        this->particle_projector.fluid_force_on_particles_two_way_coupling
          .particle_field_solution,
        std::vector<std::string>(dim, "Particle_two_way_coupling_force"),
        std::vector<DataComponentInterpretation::DataComponentInterpretation>(
          dim, DataComponentInterpretation::component_is_part_of_vector));
#endif
    }


  return solution_output_structs;
}

template <int dim>
void
FluidDynamicsVANS<dim>::monitor_mass_conservation()
{
  QGauss<dim> quadrature_formula(this->number_quadrature_points);

  FEValues<dim> fe_values(*this->mapping,
                          *this->fe,
                          quadrature_formula,
                          update_values | update_quadrature_points |
                            update_JxW_values | update_gradients |
                            update_hessians);

  FEValues<dim> fe_values_void_fraction(*this->mapping,
                                        *this->particle_projector.fe,
                                        quadrature_formula,
                                        update_values |
                                          update_quadrature_points |
                                          update_JxW_values | update_gradients);

  const FEValuesExtractors::Vector velocities(0);

  const unsigned int n_q_points = quadrature_formula.size();

  std::vector<double>         present_void_fraction_values(n_q_points);
  std::vector<Tensor<1, dim>> present_void_fraction_gradients(n_q_points);
  // Values at previous time step for transient schemes for void
  // fraction
  std::vector<double> p1_void_fraction_values(n_q_points);
  std::vector<double> p2_void_fraction_values(n_q_points);
  std::vector<double> p3_void_fraction_values(n_q_points);

  std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
  std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);

  // Values of the force function, the last component of which is the external
  // mass source
  std::vector<Vector<double>> rhs_force(n_q_points, Vector<double>(dim + 1));

  double continuity           = 0;
  double max_local_continuity = 0;
  double local_mass_source    = 0;

  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  const auto scheme = this->simulation_control->get_assembly_method();
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();

  for (const auto &cell : this->dof_handler->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          typename DoFHandler<dim>::active_cell_iterator void_fraction_cell(
            &(*this->triangulation),
            cell->level(),
            cell->index(),
            &this->particle_projector.dof_handler);
          fe_values_void_fraction.reinit(void_fraction_cell);

          // Gather void fraction (values, gradient)
          fe_values_void_fraction.get_function_values(
            particle_projector.void_fraction_locally_relevant,
            present_void_fraction_values);
          fe_values_void_fraction.get_function_gradients(
            particle_projector.void_fraction_locally_relevant,
            present_void_fraction_gradients);

          fe_values.reinit(cell);

          // Gather velocity (values and gradient)
          auto &evaluation_point = this->evaluation_point;
          fe_values[velocities].get_function_values(evaluation_point,
                                                    present_velocity_values);
          fe_values[velocities].get_function_gradients(
            evaluation_point, present_velocity_gradients);

          // Gather the external mass source at the quadrature point location
          std::vector<Point<dim>> quadrature_points =
            fe_values.get_quadrature_points();

          this->forcing_function->vector_value_list(quadrature_points,
                                                    rhs_force);

          // Gather the previous time steps depending on the number of
          // stages of the time integration scheme for the void fraction

          if (scheme !=
              Parameters::SimulationControl::TimeSteppingMethod::steady)
            {
              fe_values_void_fraction.get_function_values(
                particle_projector.previous_void_fraction[0],
                p1_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                fe_values_void_fraction.get_function_values(
                  particle_projector.previous_void_fraction[1],
                  p2_void_fraction_values);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                fe_values_void_fraction.get_function_values(
                  particle_projector.previous_void_fraction[2],
                  p3_void_fraction_values);
            }

          local_mass_source = 0;

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Calculate external mass source
              const unsigned int component_mass =
                this->fe->system_to_component_index(dim).first;
              double external_mass_source = rhs_force[q](component_mass);

              // Calculate the divergence of the velocity
              const double present_velocity_divergence =
                trace(present_velocity_gradients[q]);

              // Evaluation of global mass conservation
              local_mass_source =
                (present_velocity_values[q] *
                   present_void_fraction_gradients[q] +
                 present_void_fraction_values[q] * present_velocity_divergence -
                 external_mass_source) *
                fe_values_void_fraction.JxW(q);

              if (scheme ==
                    Parameters::SimulationControl::TimeSteppingMethod::bdf1 ||
                  scheme == Parameters::SimulationControl::TimeSteppingMethod::
                              steady_bdf)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf2)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q] +
                   bdf_coefs[2] * p2_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              if (scheme ==
                  Parameters::SimulationControl::TimeSteppingMethod::bdf3)
                local_mass_source +=
                  (bdf_coefs[0] * present_void_fraction_values[q] +
                   bdf_coefs[1] * p1_void_fraction_values[q] +
                   bdf_coefs[2] * p2_void_fraction_values[q] +
                   bdf_coefs[3] * p3_void_fraction_values[q]) *
                  fe_values_void_fraction.JxW(q);

              continuity += local_mass_source;
            }

          max_local_continuity =
            std::max(max_local_continuity, abs(local_mass_source));
        }
    }

  continuity = Utilities::MPI::sum(continuity, this->mpi_communicator);
  max_local_continuity =
    Utilities::MPI::max(max_local_continuity, this->mpi_communicator);

  this->pcout << std::setprecision(
    this->simulation_control->get_log_precision());

  this->pcout << "Global continuity equation error: " << continuity << " s^-1"
              << std::endl;
  this->pcout << "Max local continuity error: " << max_local_continuity
              << " s^-1" << std::endl;
}

template <int dim>
void
FluidDynamicsVANS<dim>::solve()
{
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

  this->setup_dofs();

  this->set_initial_condition(
    this->cfd_dem_simulation_parameters.cfd_parameters.initial_condition->type,
    this->cfd_dem_simulation_parameters.cfd_parameters.restart_parameters
      .restart);

  particle_handler.exchange_ghost_particles(true);

  while (this->simulation_control->integrate())
    {
      this->simulation_control->print_progression(this->pcout);

      // We allow the physics to update their boundary conditions
      // according to their own parameters
      this->update_boundary_conditions();
      this->multiphysics->update_boundary_conditions();

      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          vertices_cell_mapping();
          particle_projector.initialize_void_fraction(
            this->simulation_control->get_current_time());
          this->iterate();
        }
      else
        {
          NavierStokesBase<dim, GlobalVectorType, IndexSet>::refine_mesh();
          vertices_cell_mapping();
          calculate_void_fraction(this->simulation_control->get_current_time());
          this->iterate();
        }

      this->postprocess(false);
      monitor_mass_conservation();
      finish_time_step_fd();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D solver to ensure that the
// library is valid before we actually compile the solver
template class FluidDynamicsVANS<2>;
template class FluidDynamicsVANS<3>;
