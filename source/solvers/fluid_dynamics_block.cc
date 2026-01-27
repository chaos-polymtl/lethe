// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/fluid_dynamics_block.h>
#include <solvers/isothermal_compressible_navier_stokes_vof_assembler.h>
#include <solvers/navier_stokes_assemblers.h>
#include <solvers/navier_stokes_vof_assemblers.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/numerics/vector_tools.h>

// Constructor for class FluidDynamicsBlock
template <int dim>
FluidDynamicsBlock<dim>::FluidDynamicsBlock(
  SimulationParameters<dim> &p_nsparam)
  : NavierStokesBase<dim, GlobalBlockVectorType, std::vector<IndexSet>>(
      p_nsparam)
{}

template <int dim>
FluidDynamicsBlock<dim>::~FluidDynamicsBlock()
{
  this->dof_handler->clear();
}

template <int dim>
void
FluidDynamicsBlock<dim>::setup_assemblers()
{
  this->assemblers.clear();

  // Buoyant force
  if (this->simulation_parameters.multiphysics.buoyancy_force)
    {
      this->assemblers.emplace_back(std::make_shared<BuoyancyAssembly<dim>>(
        this->simulation_control,
        this->simulation_parameters.physical_properties_manager
          .get_reference_temperature()));
    }

  // ALE
  if (this->simulation_parameters.ale.enabled())
    {
      this->assemblers.emplace_back(
        std::make_shared<NavierStokesAssemblerALE<dim>>(
          this->simulation_control, this->simulation_parameters.ale));
    }

  if (this->simulation_parameters.multiphysics.VOF)
    {
      // Time-stepping schemes
      if (time_stepping_is_bdf(
            this->simulation_control->get_assembly_method()) &&
          this->simulation_parameters.physical_properties_manager
            .density_is_constant())
        {
          this->assemblers.emplace_back(
            std::make_shared<GLSNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }
      else if (time_stepping_is_bdf(
                 this->simulation_control->get_assembly_method()))
        {
          this->assemblers.emplace_back(
            std::make_shared<
              GLSIsothermalCompressibleNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }

      // Darcy force for phase change simulations
      if (this->simulation_parameters.velocity_sources.darcy_type ==
          Parameters::VelocitySource::DarcySourceType::phase_change)
        {
          AssertThrow(this->simulation_parameters.multiphysics.heat_transfer,
                      PhaseChangeDarcyModelRequiresTemperature());
          this->assemblers.emplace_back(
            std::make_shared<PhaseChangeDarcyVOFAssembler<dim>>(
              this->simulation_parameters.physical_properties_manager
                .get_phase_change_parameters_vector()));
        }

      if (!this->simulation_parameters.physical_properties_manager
             .density_is_constant())
        {
          this->assemblers.emplace_back(
            std::make_shared<
              GLSIsothermalCompressibleNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
      else
        {
          // Core assembler
          this->assemblers.emplace_back(
            std::make_shared<GLSNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
    }
  else
    {
      // Time-stepping schemes
      if (time_stepping_is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.emplace_back(
            std::make_shared<GLSNavierStokesAssemblerBDF<dim>>(
              this->simulation_control));
        }

      // Velocity sources term
      if (this->simulation_parameters.velocity_sources.rotating_frame_type ==
          Parameters::VelocitySource::RotatingFrameType::srf)
        {
          this->assemblers.emplace_back(
            std::make_shared<GLSNavierStokesAssemblerSRF<dim>>(
              this->simulation_parameters.velocity_sources));
        }

      // Darcy force for phase change simulations in single phase
      if (this->simulation_parameters.velocity_sources.darcy_type ==
          Parameters::VelocitySource::DarcySourceType::phase_change)
        {
          AssertThrow(this->simulation_parameters.multiphysics.heat_transfer,
                      PhaseChangeDarcyModelRequiresTemperature());
          this->assemblers.emplace_back(
            std::make_shared<PhaseChangeDarcyAssembly<dim>>(
              this->simulation_parameters.physical_properties_manager
                .get_physical_properties_parameters()
                .fluids[0]
                .phase_change_parameters));
        }

      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with non-Newtonian viscosity
          this->assemblers.emplace_back(
            std::make_shared<BlockNavierStokesAssemblerNonNewtonianCore<dim>>(
              this->simulation_control, gamma));
        }
      else
        {
          // Core default assembler
          if ((this->simulation_parameters.stabilization
                 .use_default_stabilization == true) ||
              this->simulation_parameters.stabilization.stabilization ==
                Parameters::Stabilization::NavierStokesStabilization::grad_div)
            this->assemblers.emplace_back(
              std::make_shared<BlockNavierStokesAssemblerCore<dim>>(
                this->simulation_control, gamma));

          else
            throw std::runtime_error(
              "Using the GD solver with a stabilization other than the grad_div "
              "stabilization will lead to an unstable block solver that is unable to converge");
        }
    }
}

template <int dim>
void
FluidDynamicsBlock<dim>::update_multiphysics_time_average_solution()
{
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->multiphysics->set_block_time_average_solution(
        PhysicsID::fluid_dynamics,
        this->average_velocities->get_average_velocities());
    }
}

template <int dim>
void
FluidDynamicsBlock<dim>::assemble_system_matrix()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");
  // this->simulation_control->set_assembly_method(this->time_stepping_method);

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
    }


  WorkStream::run(
    this->dof_handler->begin_active(),
    this->dof_handler->end(),
    *this,
    &FluidDynamicsBlock::assemble_local_system_matrix,
    &FluidDynamicsBlock::copy_local_matrix_to_global_matrix,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));


  system_matrix.compress(VectorOperation::add);

  // Finally we move pressure mass matrix into a separate matrix:
  pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));
  pressure_mass_matrix.copy_from(system_matrix.block(1, 1));

  // Note that settings this pressure block to zero is not identical to
  // not assembling anything in this block, because this operation here
  // will (incorrectly) delete diagonal entries that come in from
  // hanging node constraints for pressure DoFs. This means that our
  // whole system matrix will have rows that are completely
  // zero. Luckily, FGMRES handles these rows without any problem.
  system_matrix.block(1, 1) = 0;
}



template <int dim>
void
FluidDynamicsBlock<dim>::assemble_local_system_matrix(
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
  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> &dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        this->multiphysics->get_solution(PhysicsID::VOF),
        this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> &dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht.get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();


  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_matrix(scratch_data, copy_data);
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
FluidDynamicsBlock<dim>::copy_local_matrix_to_global_matrix(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_matrix,
                                              copy_data.local_dof_indices,
                                              system_matrix);
}


template <int dim>
void
FluidDynamicsBlock<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");
  // this->simulation_control->set_assembly_method(this->time_stepping_method);

  this->system_rhs = 0;
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
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> &dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht.get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }


  WorkStream::run(
    this->dof_handler->begin_active(),
    this->dof_handler->end(),
    *this,
    &FluidDynamicsBlock::assemble_local_system_rhs,
    &FluidDynamicsBlock::copy_local_rhs_to_global_rhs,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}

template <int dim>
void
FluidDynamicsBlock<dim>::setup_preconditioner()
{
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .preconditioner == Parameters::LinearSolver::PreconditionerType::ilu)
    {
      setup_ILU();
    }
  else if (this->simulation_parameters.linear_solver
             .at(PhysicsID::fluid_dynamics)
             .preconditioner ==
           Parameters::LinearSolver::PreconditionerType::amg)
    {
      setup_AMG();
    }
  else
    AssertThrow(
      false,
      ExcMessage(
        "This linear solver does not support this preconditioner. Only <ilu> and <amg> preconditioners are supported."));
}

template <int dim>
void
FluidDynamicsBlock<dim>::assemble_local_system_rhs(
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

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> &dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        this->multiphysics->get_solution(PhysicsID::VOF),
        this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        this->multiphysics->get_previous_solutions(PhysicsID::VOF));
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> &dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      typename DoFHandler<dim>::active_cell_iterator temperature_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        &dof_handler_ht);

      scratch_data.reinit_heat_transfer(
        temperature_cell,
        this->multiphysics->get_solution(PhysicsID::heat_transfer),
        this->multiphysics->get_previous_solutions(PhysicsID::heat_transfer));
    }

  scratch_data.calculate_physical_properties();
  copy_data.reset();
  for (auto &assembler : this->assemblers)
    {
      assembler->assemble_rhs(scratch_data, copy_data);
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}



template <int dim>
void
FluidDynamicsBlock<dim>::copy_local_rhs_to_global_rhs(
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
void
FluidDynamicsBlock<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  auto &system_rhs = this->system_rhs;
  system_rhs       = 0;
  FEValues<dim>               fe_values(*this->mapping,
                          *this->fe,
                          *this->cell_quadrature,
                          update_values | update_quadrature_points |
                            update_JxW_values);
  const unsigned int          dofs_per_cell = this->fe->dofs_per_cell;
  const unsigned int          n_q_points    = this->cell_quadrature->size();
  FullMatrix<double>          local_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>              local_rhs(dofs_per_cell);
  std::vector<Vector<double>> initial_velocity(n_q_points,
                                               Vector<double>(dim + 1));
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  const FEValuesExtractors::Vector     velocities(0);
  const FEValuesExtractors::Scalar     pressure(dim);

  Tensor<1, dim> rhs_initial_velocity_pressure;
  double         rhs_initial_pressure;

  std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
  std::vector<double>         phi_p(dofs_per_cell);

  for (const auto &cell : this->dof_handler->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          local_matrix = 0;
          local_rhs    = 0;
          this->simulation_parameters.initial_condition->uvwp.vector_value_list(
            fe_values.get_quadrature_points(), initial_velocity);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  phi_p[k] = fe_values[pressure].value(k, q);
                  phi_u[k] = fe_values[velocities].value(k, q);
                }

              // Establish the rhs tensor operator
              for (int i = 0; i < dim; ++i)
                {
                  const unsigned int component_i =
                    this->fe->system_to_component_index(i).first;
                  rhs_initial_velocity_pressure[i] =
                    initial_velocity[q](component_i);
                }
              rhs_initial_pressure = initial_velocity[q](dim);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Matrix assembly
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi_u[j] * phi_u[i]) * fe_values.JxW(q);
                      local_matrix(i, j) +=
                        (phi_p[j] * phi_p[i]) * fe_values.JxW(q);
                    }
                  local_rhs(i) += (phi_u[i] * rhs_initial_velocity_pressure +
                                   phi_p[i] * rhs_initial_pressure) *
                                  fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          this->nonzero_constraints.distribute_local_to_global(
            local_matrix,
            local_rhs,
            local_dof_indices,
            system_matrix,
            system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}


template <int dim>
void
FluidDynamicsBlock<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "Setup DOFs");

  system_matrix.clear();

  this->dof_handler->distribute_dofs(*this->fe);
  // DoFRenumbering::Cuthill_McKee(this->dof_handler);


  std::vector<unsigned int> block_component(dim + 1, 0);
  block_component[dim] = 1;
  DoFRenumbering::component_wise(*this->dof_handler, block_component);


  // To be used to replace the above part once 9.2 release is out and TRAVIS-CI
  // version is updated
  dofs_per_block =
    DoFTools::count_dofs_per_fe_block(*this->dof_handler, block_component);

  unsigned int dof_u = dofs_per_block[0];
  unsigned int dof_p = dofs_per_block[1];

  this->locally_owned_dofs.resize(2);
  this->locally_owned_dofs[0] =
    this->dof_handler->locally_owned_dofs().get_view(0, dof_u);
  this->locally_owned_dofs[1] =
    this->dof_handler->locally_owned_dofs().get_view(dof_u, dof_u + dof_p);

  IndexSet locally_relevant_dofs_acquisition;
  locally_relevant_dofs_acquisition =
    DoFTools::extract_locally_relevant_dofs(*this->dof_handler);
  this->locally_relevant_dofs.resize(2);
  this->locally_relevant_dofs[0] =
    locally_relevant_dofs_acquisition.get_view(0, dof_u);
  this->locally_relevant_dofs[1] =
    locally_relevant_dofs_acquisition.get_view(dof_u, dof_u + dof_p);

  // Non-zero constraints
  this->define_non_zero_constraints();

  // Check whether the boundary conditions specified in the parameter file are
  // available for this solver
  for (auto const &[id, type] :
       this->simulation_parameters.boundary_conditions.type)
    {
      if (type == BoundaryConditions::BoundaryType::pressure ||
          type == BoundaryConditions::BoundaryType::function_weak ||
          type == BoundaryConditions::BoundaryType::partial_slip)
        {
          Assert(
            false,
            ExcMessage(
              "The following boundary conditions are not supported by the lethe-fluid-block application: pressure, function weak, partial slip and outlet."));
        }
    }

  // Zero constraints
  this->define_zero_constraints();

  // Operations on the following vectors (addition,
  // multiplication, etc.) can only be done if these are reinitialized WITHOUT
  // locally_relevant_dofs (i.e. without ghost DoFs). This is why most of them
  // are reinitialized both with and without locally_relevant_dofs.
  this->present_solution->reinit(this->locally_owned_dofs,
                                 this->locally_relevant_dofs,
                                 this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);
  this->evaluation_point.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  // Initialize vector of previous solutions
  for (auto &solution : *this->previous_solutions)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      this->mpi_communicator);
    }


  if (this->simulation_control->is_sdirk())
    {
      // Reinitialize vectors used for the SDIRK methods
      this->sdirk_vectors.sum_bi_ki.reinit(this->locally_owned_dofs,
                                           this->locally_relevant_dofs,
                                           this->mpi_communicator);
      this->sdirk_vectors.local_sum_bi_ki.reinit(this->locally_owned_dofs,
                                                 this->mpi_communicator);
      this->sdirk_vectors.sum_over_previous_stages.reinit(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->mpi_communicator);
      this->sdirk_vectors.local_sum_over_previous_stages.reinit(
        this->locally_owned_dofs, this->mpi_communicator);
      this->sdirk_vectors.locally_owned_for_calculation.reinit(
        this->locally_owned_dofs, this->mpi_communicator);

      for (auto &solution : this->sdirk_vectors.previous_k_j_solutions)
        {
          solution.reinit(this->locally_owned_dofs,
                          this->locally_relevant_dofs,
                          this->mpi_communicator);
        }
    }


  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);


  sparsity_pattern.reinit(this->locally_owned_dofs,
                          this->locally_owned_dofs,
                          this->locally_relevant_dofs,
                          MPI_COMM_WORLD);

  Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
  for (int c = 0; c < dim + 1; ++c)
    for (int d = 0; d < dim + 1; ++d)
      if (!((c == dim) && (d == dim)))
        coupling[c][d] = DoFTools::always;
      else
        coupling[c][d] = DoFTools::always;

  DoFTools::make_sparsity_pattern(*this->dof_handler,
                                  coupling,
                                  sparsity_pattern,
                                  this->nonzero_constraints,
                                  true,
                                  Utilities::MPI::this_mpi_process(
                                    MPI_COMM_WORLD));

  sparsity_pattern.compress();

  system_matrix.reinit(sparsity_pattern);
  pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1));

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      AssertThrow(this->simulation_parameters.mesh_adaptation.type ==
                    Parameters::MeshAdaptation::Type::none,
                  ExcMessage(
                    "Time-averaging velocities and calculating reynolds "
                    "stresses are currently unavailable for mesh "
                    "adaptation."));

      this->average_velocities->initialize_vectors(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->fe->n_dofs_per_vertex(),
        this->mpi_communicator);
    }


  double global_volume =
    GridTools::volume(*this->triangulation, *this->mapping);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler->n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;

  this->multiphysics->set_dof_handler(PhysicsID::fluid_dynamics,
                                      this->dof_handler);
  this->multiphysics->set_block_solution(PhysicsID::fluid_dynamics,
                                         this->present_solution);
  this->multiphysics->set_block_previous_solutions(PhysicsID::fluid_dynamics,
                                                   this->previous_solutions);
}


template <int dim>
void
FluidDynamicsBlock<dim>::set_solution_vector(double value)
{
  *this->present_solution = value;
}


/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
FluidDynamicsBlock<dim>::set_initial_condition_fd(
  Parameters::FluidDynamicsInitialConditionType initial_condition_type,
  bool                                          restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::L2projection)
    {
      assemble_L2_projection();
      solve_L2_system(1e-15, 1e-15);
      *this->present_solution = this->newton_update;
      this->finish_time_step();
    }
  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->finish_time_step();
    }

  else if (initial_condition_type ==
           Parameters::FluidDynamicsInitialConditionType::viscous)
    {
      this->set_nodal_values();
      std::shared_ptr<RheologicalModel> original_viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Temporarily set the rheology to be newtonian with predefined viscosity
      std::shared_ptr<Newtonian> temporary_rheology =
        std::make_shared<Newtonian>(
          this->simulation_parameters.initial_condition->kinematic_viscosity);

      this->simulation_parameters.physical_properties_manager.set_rheology(
        temporary_rheology);


      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::steady);
      PhysicsSolver<GlobalBlockVectorType>::solve_governing_system();
      this->finish_time_step();

      this->simulation_parameters.physical_properties_manager.set_rheology(
        original_viscosity_model);
    }
  else
    {
      throw std::runtime_error("GDNS - Initial condition could not be set");
    }
}

template <int dim>
void
FluidDynamicsBlock<dim>::solve_linear_system()
{
  const double absolute_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .relative_residual;

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .solver == Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(absolute_residual, relative_residual);
  else
    AssertThrow(
      this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .solver == Parameters::LinearSolver::SolverType::gmres,
      ExcMessage(
        "This linear solver is not allowed. Only <gmres> linear solver is allowed."));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
FluidDynamicsBlock<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const unsigned int ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;

  velocity_ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();
  pressure_ilu_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionILU>();

  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  velocity_ilu_preconditioner->initialize(system_matrix.block(0, 0),
                                          preconditionerOptions);

  pressure_ilu_preconditioner->initialize(system_matrix.block(1, 1),
                                          preconditionerOptions);
  system_ilu_preconditioner = std::make_shared<
    BlockSchurPreconditioner<TrilinosWrappers::PreconditionILU>>(
    gamma,
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_ilu_preconditioner),
    &(*pressure_ilu_preconditioner),
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics));
}

template <int dim>
void
FluidDynamicsBlock<dim>::setup_AMG()
{
  TimerOutput::Scope t(this->computing_timer, "Setup AMG");

  //**********************************************
  // Trillinos Wrapper AMG Preconditioner
  //*********************************************

  // Constant modes for velocity
  std::vector<std::vector<bool>> velocity_constant_modes;
  FEValuesExtractors::Vector     velocities(0);

  ComponentMask velocity_components = this->fe->component_mask(velocities);
  velocity_constant_modes =
    DoFTools::extract_constant_modes(*this->dof_handler, velocity_components);

  // Constant modes for pressure
  std::vector<std::vector<bool>> pressure_constant_modes;
  FEValuesExtractors::Scalar     pressure(dim);
  ComponentMask pressure_components = this->fe->component_mask(pressure);
  pressure_constant_modes =
    DoFTools::extract_constant_modes(*this->dof_handler, pressure_components);

  this->computing_timer.enter_subsection("AMG_velocity");
  const bool elliptic_velocity     = false;
  bool       higher_order_elements = false;
  if (this->fe->degree > 1)
    higher_order_elements = true;
  const unsigned int n_cycles =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_n_cycles;
  const bool w_cycle =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_w_cycles;
  const double aggregation_threshold =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .amg_smoother_overlap;
  const bool  output_details = false;
  const char *smoother_type  = "Chebyshev";  //"ILU";
  const char *coarse_type    = "Amesos-KLU"; //"ILU";

  velocity_amg_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionAMG>();
  pressure_amg_preconditioner =
    std::make_shared<TrilinosWrappers::PreconditionAMG>();

  TrilinosWrappers::PreconditionAMG::AdditionalData
    velocity_preconditioner_options(elliptic_velocity,
                                    higher_order_elements,
                                    n_cycles,
                                    w_cycle,
                                    aggregation_threshold,
                                    velocity_constant_modes,
                                    smoother_sweeps,
                                    smoother_overlap,
                                    output_details,
                                    smoother_type,
                                    coarse_type);

  Teuchos::ParameterList              velocity_parameter_ml;
  std::unique_ptr<Epetra_MultiVector> velocity_distributed_constant_modes;
  velocity_preconditioner_options.set_parameters(
    velocity_parameter_ml,
    velocity_distributed_constant_modes,
    system_matrix.block(0, 0));
  velocity_amg_preconditioner->initialize(system_matrix.block(0, 0),
                                          velocity_parameter_ml);
  this->computing_timer.leave_subsection("AMG_velocity");

  this->computing_timer.enter_subsection("AMG_pressure");
  const bool elliptic_pressure = true;
  higher_order_elements        = false;
  if (this->pressure_fem_degree > 1)
    higher_order_elements = true;
  TrilinosWrappers::PreconditionAMG::AdditionalData
                         pressure_preconditioner_options(elliptic_pressure,
                                    higher_order_elements,
                                    n_cycles,
                                    w_cycle,
                                    aggregation_threshold,
                                    pressure_constant_modes,
                                    smoother_sweeps,
                                    smoother_overlap,
                                    output_details,
                                    smoother_type,
                                    coarse_type);
  Teuchos::ParameterList pressure_parameter_ml;
  std::unique_ptr<Epetra_MultiVector> pressure_distributed_constant_modes;
  velocity_preconditioner_options.set_parameters(
    pressure_parameter_ml,
    pressure_distributed_constant_modes,
    system_matrix.block(0, 0));
  pressure_amg_preconditioner->initialize(system_matrix.block(1, 1),
                                          pressure_parameter_ml);
  this->computing_timer.leave_subsection("AMG_pressure");


  GlobalBlockVectorType completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);


  system_amg_preconditioner = std::make_shared<
    BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG>>(
    gamma,
    this->simulation_parameters.physical_properties_manager
      .get_kinematic_viscosity_scale(),
    system_matrix,
    pressure_mass_matrix,
    &(*velocity_amg_preconditioner),
    &(*pressure_amg_preconditioner),
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics));
}



template <int dim>
void
FluidDynamicsBlock<dim>::solve_system_GMRES(const double absolute_residual,
                                            const double relative_residual)
{
  const AffineConstraints<double> &zero_constraints = this->zero_constraints;
  const double rescale_metric   = this->get_residual_rescale_metric();
  const double current_residual = this->system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * current_residual, absolute_residual);
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_tolerance * rescale_metric;

  GlobalBlockVectorType completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               non_rescaled_linear_solver_tolerance,
                               true,
                               true);

  SolverFGMRES<GlobalBlockVectorType> solver(solver_control);

  setup_preconditioner();

  {
    TimerOutput::Scope t(this->computing_timer, "Solve linear system");
    if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .preconditioner == Parameters::LinearSolver::PreconditionerType::ilu)
      solver.solve(this->system_matrix,
                   this->newton_update,
                   this->system_rhs,
                   *system_ilu_preconditioner);
    else if (this->simulation_parameters.linear_solver
               .at(PhysicsID::fluid_dynamics)
               .preconditioner ==
             Parameters::LinearSolver::PreconditionerType::amg)
      solver.solve(this->system_matrix,
                   this->newton_update,
                   this->system_rhs,
                   *system_amg_preconditioner);
    else
      AssertThrow(
        this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::ilu ||
          this->simulation_parameters.linear_solver
              .at(PhysicsID::fluid_dynamics)
              .preconditioner ==
            Parameters::LinearSolver::PreconditionerType::amg,
        ExcMessage(
          "This linear solver does not support this preconditioner. Only <ilu> and <amg> preconditioners are supported."));

    if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
          .verbosity != Parameters::Verbosity::quiet)
      {
        this->pcout << "  -Iterative solver took : "
                    << solver_control.last_step()
                    << " steps to reach a residual norm of "
                    << solver_control.last_value() / rescale_metric
                    << std::endl;
      }

    zero_constraints.distribute(this->newton_update);
  }
}

template <int dim>
void
FluidDynamicsBlock<dim>::solve_L2_system(double absolute_residual,
                                         double relative_residual)
{
  TimerOutput::Scope t(this->computing_timer, "Solve linear system");

  const AffineConstraints<double> &nonzero_constraints =
    this->nonzero_constraints;
  const double rescale_metric   = this->get_residual_rescale_metric();
  const double current_residual = this->system_rhs.l2_norm() / rescale_metric;
  const double linear_solver_tolerance =
    std::max(relative_residual * current_residual, absolute_residual);
  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  const double non_rescaled_linear_solver_tolerance =
    linear_solver_tolerance * rescale_metric;

  GlobalBlockVectorType completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(this->simulation_parameters.linear_solver
                                 .at(PhysicsID::fluid_dynamics)
                                 .max_iterations,
                               non_rescaled_linear_solver_tolerance,
                               true,
                               true);
  SolverFGMRES<GlobalBlockVectorType> solver(solver_control);
  TrilinosWrappers::PreconditionILU   pmass_preconditioner;

  //**********************************************
  // Trillinos Wrapper ILU Preconditioner
  //*********************************************
  const unsigned int ilu_fill =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_fill;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
      .ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    ilu_fill, ilu_atol, ilu_rtol, 0);
  pmass_preconditioner.initialize(pressure_mass_matrix, preconditionerOptions);


  const BlockDiagPreconditioner<TrilinosWrappers::PreconditionILU>
    preconditioner(system_matrix, pmass_preconditioner, solver_control);

  solver.solve(system_matrix,
               completely_distributed_solution,
               this->system_rhs,
               preconditioner);

  if (this->simulation_parameters.linear_solver.at(PhysicsID::fluid_dynamics)
        .verbosity != Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
    }

  nonzero_constraints.distribute(completely_distributed_solution);
  this->newton_update = completely_distributed_solution;
}

template <int dim>
void
FluidDynamicsBlock<dim>::multi_stage_preresolution(
  unsigned int                                      stage,
  Parameters::SimulationControl::TimeSteppingMethod method)
{
  // Copy the reference to some of the vectors to enhance readability below
  auto &previous_k_j_solutions   = this->sdirk_vectors.previous_k_j_solutions;
  auto &sum_over_previous_stages = this->sdirk_vectors.sum_over_previous_stages;
  auto &local_sum_over_previous_stages =
    this->sdirk_vectors.local_sum_over_previous_stages;
  auto &locally_owned_for_calculation =
    this->sdirk_vectors.locally_owned_for_calculation;

  // If a SDIRK method is selected, we need to solve as many
  // nonlinear systems as the number of stages.

  // At each stage, the value of the (a_ij), (b_i) and (c_i)
  // coefficients are differents

  SDIRKTable     table = sdirk_table(method);
  SDIRKStageData stage_data(table, stage + 1);
  const double   a_ii = stage_data.a_ij[stage];

  // At each stage, we need to recompute sum(a_{ij} * k_j)
  // = sum_over_previous_stages
  local_sum_over_previous_stages = 0;

  if (stage > 0)
    {
      // At the first stage, the sum_over_previous_stages is set to
      // 0. But for the next stages, we need to update this
      // sum_over_previous_stages
      for (unsigned int p = 0; p < stage; ++p)
        {
          locally_owned_for_calculation = previous_k_j_solutions[p];
          local_sum_over_previous_stages.add(stage_data.a_ij[p] / a_ii,
                                             locally_owned_for_calculation);
        }
    }

  sum_over_previous_stages = local_sum_over_previous_stages;
}

template <int dim>
void
FluidDynamicsBlock<dim>::multi_stage_postresolution(
  unsigned int                                      stage,
  Parameters::SimulationControl::TimeSteppingMethod method,
  double                                            time_step)
{
  // Copy the reference to some of the vectors to enhance readability below
  auto &present_solution       = *this->present_solution;
  auto &previous_solutions     = *this->previous_solutions;
  auto &local_evaluation_point = this->local_evaluation_point;
  auto &local_sum_over_previous_stages =
    this->sdirk_vectors.local_sum_over_previous_stages;
  auto &locally_owned_for_calculation =
    this->sdirk_vectors.locally_owned_for_calculation;
  auto &previous_k_j_solutions = this->sdirk_vectors.previous_k_j_solutions;
  auto &sum_bi_ki              = this->sdirk_vectors.sum_bi_ki;
  auto &local_sum_bi_ki        = this->sdirk_vectors.local_sum_bi_ki;

  SDIRKTable     table = sdirk_table(method);
  SDIRKStageData stage_data(table, stage + 1);
  const double   a_ii = stage_data.a_ij[stage];

  // Once we have solved the nonlinear system for the velocity, we
  // want to store the value of the coefficient k_i for the final
  // sum b_i*k_i with k_i = (u*_{i} - u_{n})/(time_step*a_ii) -
  // sum_over_previous_stages
  locally_owned_for_calculation = previous_solutions[0];
  local_evaluation_point        = present_solution;
  local_evaluation_point.add(-1.0, locally_owned_for_calculation);
  local_evaluation_point *= 1.0 / (time_step * a_ii);
  local_evaluation_point.add(-1.0, local_sum_over_previous_stages);

  // We store the value of the present_k_i_solution in the
  // previous_k_j_solutions vector.
  previous_k_j_solutions[stage] = local_evaluation_point;

  // We update the sum of b_i*k_i
  const double b_i = stage_data.b_i;
  local_sum_bi_ki.add(b_i, local_evaluation_point);
  sum_bi_ki = local_sum_bi_ki;
}

template <int dim>
void
FluidDynamicsBlock<dim>::update_multi_stage_solution(double time_step)
{
  // Copy the reference to some of the vectors to enhance readability below
  auto &present_solution       = *this->present_solution;
  auto &previous_solutions     = *this->previous_solutions;
  auto &local_evaluation_point = this->local_evaluation_point;
  auto &sum_bi_ki              = this->sdirk_vectors.sum_bi_ki;
  auto &local_sum_bi_ki        = this->sdirk_vectors.local_sum_bi_ki;

  // At each time iteration, we update the value of present_solution
  local_sum_bi_ki        = sum_bi_ki;
  local_evaluation_point = previous_solutions[0];
  local_evaluation_point.add(time_step, local_sum_bi_ki);
  present_solution = local_evaluation_point;
}

/*
 * Generic CFD Solver application
 * Handles the majority of the cases for the GD-NS solver
 */
template <int dim>
void
FluidDynamicsBlock<dim>::solve()
{
  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();
  this->box_refine_mesh(this->simulation_parameters.restart_parameters.restart);
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);
  this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          for (unsigned int i = 0;
               i <
               this->simulation_parameters.mesh_adaptation.initial_refinement;
               i++)
            NavierStokesBase<dim,
                             GlobalBlockVectorType,
                             std::vector<IndexSet>>::refine_mesh();

          this->iterate();
        }
      else
        {
          NavierStokesBase<dim, GlobalBlockVectorType, std::vector<IndexSet>>::
            refine_mesh();
          this->iterate();
        }
      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class FluidDynamicsBlock<2>;
template class FluidDynamicsBlock<3>;
