/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Bruno Blais, Polytechnique Montreal, 2019-
 */

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/multiphysics.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/gls_navier_stokes.h>
#include <solvers/isothermal_compressible_navier_stokes_assembler.h>
#include <solvers/navier_stokes_vof_assemblers.h>

#include <deal.II/base/work_stream.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/numerics/vector_tools.h>



// Constructor for class GLSNavierStokesSolver
template <int dim>
GLSNavierStokesSolver<dim>::GLSNavierStokesSolver(
  SimulationParameters<dim> &p_nsparam)
  : NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>(p_nsparam)
{
  initial_preconditioner_fill_level =
    ((this->simulation_parameters.linear_solver.solver ==
      Parameters::LinearSolver::SolverType::amg) ?
       this->simulation_parameters.linear_solver.amg_precond_ilu_fill :
       this->simulation_parameters.linear_solver.ilu_precond_fill);
}

template <int dim>
GLSNavierStokesSolver<dim>::~GLSNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_dofs_fd()
{
  TimerOutput::Scope t(this->computing_timer, "setup_dofs");

  // Clear the preconditioner before the matrix they are associated with is
  // cleared
  amg_preconditioner.reset();
  ilu_preconditioner.reset();
  current_preconditioner_fill_level = initial_preconditioner_fill_level;


  // Now reset system matrix
  system_matrix.clear();

  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

  FEValuesExtractors::Vector velocities(0);

  // Non Zero constraints
  define_non_zero_constraints();

  // Zero constraints
  define_zero_constraints();

  this->present_solution.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  this->evaluation_point.reinit(this->locally_owned_dofs,
                                this->locally_relevant_dofs,
                                this->mpi_communicator);

  // Initialize vector of previous solutions
  for (auto &solution : this->previous_solutions)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      this->mpi_communicator);
    }

  // If SDIRK type of methods are used, initialize solution stages
  for (auto &solution : this->solution_stages)
    {
      solution.reinit(this->locally_owned_dofs,
                      this->locally_relevant_dofs,
                      this->mpi_communicator);
    }

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);
  auto &                 nonzero_constraints = this->get_nonzero_constraints();
  DynamicSparsityPattern dsp(this->locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  false);
  SparsityTools::distribute_sparsity_pattern(
    dsp,
    this->dof_handler.locally_owned_dofs(),
    this->mpi_communicator,
    this->locally_relevant_dofs);

  system_matrix.reinit(this->locally_owned_dofs,
                       this->locally_owned_dofs,
                       dsp,
                       this->mpi_communicator);

  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->average_velocities->initialize_vectors(
        this->locally_owned_dofs,
        this->locally_relevant_dofs,
        this->fe->n_dofs_per_vertex(),
        this->mpi_communicator);

      if (this->simulation_parameters.restart_parameters.checkpoint)
        {
          this->average_velocities->initialize_checkpoint_vectors(
            this->locally_owned_dofs,
            this->locally_relevant_dofs,
            this->mpi_communicator);
        }
    }

  double global_volume =
    GridTools::volume(*this->triangulation, *this->mapping);

  this->pcout << "   Number of active cells:       "
              << this->triangulation->n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: "
              << this->dof_handler.n_dofs() << std::endl;
  this->pcout << "   Volume of triangulation:      " << global_volume
              << std::endl;


  // Provide the fluid dynamics dof_handler and present solution to the
  // multiphysics interface
  this->multiphysics->set_dof_handler(PhysicsID::fluid_dynamics,
                                      &this->dof_handler);
  this->multiphysics->set_solution(PhysicsID::fluid_dynamics,
                                   &this->present_solution);
  this->multiphysics->set_previous_solutions(PhysicsID::fluid_dynamics,
                                             &this->previous_solutions);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::update_multiphysics_time_average_solution()
{
  if (this->simulation_parameters.post_processing.calculate_average_velocities)
    {
      this->multiphysics->set_time_average_solution(
        PhysicsID::fluid_dynamics,
        &this->average_velocities->get_average_velocities());
    }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::update_boundary_conditions()
{
  double time = this->simulation_control->get_current_time();
  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .u.set_time(time);
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .v.set_time(time);
      this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
        .w.set_time(time);
      this->simulation_parameters.boundary_conditions.bcPressureFunction[i_bc]
        .p.set_time(time);
    }
  define_non_zero_constraints();
  // Distribute constraints
  auto &nonzero_constraints = this->nonzero_constraints;
  nonzero_constraints.distribute(this->local_evaluation_point);
  this->present_solution = this->local_evaluation_point;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::define_non_zero_constraints()
{
  double time = this->simulation_control->get_current_time();
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  // Non-zero constraints
  auto &nonzero_constraints = this->get_nonzero_constraints();
  {
    nonzero_constraints.clear();
    nonzero_constraints.reinit(this->locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(this->dof_handler,
                                            nonzero_constraints);
    for (unsigned int i_bc = 0;
         i_bc < this->simulation_parameters.boundary_conditions.size;
         ++i_bc)
      {
        if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
            BoundaryConditions::BoundaryType::noslip)
          {
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              dealii::Functions::ZeroFunction<dim>(dim + 1),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::slip)
          {
            std::set<types::boundary_id> no_normal_flux_boundaries;
            no_normal_flux_boundaries.insert(
              this->simulation_parameters.boundary_conditions.id[i_bc]);
            VectorTools::compute_no_normal_flux_constraints(
              this->dof_handler,
              0,
              no_normal_flux_boundaries,
              nonzero_constraints,
              *this->mapping);
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::function)
          {
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .u.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .v.set_time(time);
            this->simulation_parameters.boundary_conditions.bcFunctions[i_bc]
              .w.set_time(time);
            VectorTools::interpolate_boundary_values(
              *this->mapping,
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              NavierStokesFunctionDefined<dim>(
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .u,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .v,
                &this->simulation_parameters.boundary_conditions
                   .bcFunctions[i_bc]
                   .w),
              nonzero_constraints,
              this->fe->component_mask(velocities));
          }
        else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
                 BoundaryConditions::BoundaryType::periodic)
          {
            DoFTools::make_periodicity_constraints(
              this->dof_handler,
              this->simulation_parameters.boundary_conditions.id[i_bc],
              this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
              this->simulation_parameters.boundary_conditions
                .periodic_direction[i_bc],
              nonzero_constraints);
          }
      }
  }

  this->establish_solid_domain(true);

  nonzero_constraints.close();
}

template <int dim>
void
GLSNavierStokesSolver<dim>::define_zero_constraints()
{
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure(dim);
  this->zero_constraints.clear();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);
  this->zero_constraints.reinit(this->locally_relevant_dofs);

  DoFTools::make_hanging_node_constraints(this->dof_handler,
                                          this->zero_constraints);

  for (unsigned int i_bc = 0;
       i_bc < this->simulation_parameters.boundary_conditions.size;
       ++i_bc)
    {
      if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
          BoundaryConditions::BoundaryType::slip)
        {
          std::set<types::boundary_id> no_normal_flux_boundaries;
          no_normal_flux_boundaries.insert(
            this->simulation_parameters.boundary_conditions.id[i_bc]);
          VectorTools::compute_no_normal_flux_constraints(
            this->dof_handler,
            0,
            no_normal_flux_boundaries,
            this->zero_constraints,
            *this->mapping);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::periodic)
        {
          DoFTools::make_periodicity_constraints(
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            this->simulation_parameters.boundary_conditions.periodic_id[i_bc],
            this->simulation_parameters.boundary_conditions
              .periodic_direction[i_bc],
            this->zero_constraints);
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::pressure)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::function_weak)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::partial_slip)
        {
          /*do nothing*/
        }
      else if (this->simulation_parameters.boundary_conditions.type[i_bc] ==
               BoundaryConditions::BoundaryType::outlet)
        {
          /*do nothing*/
        }
      else
        {
          VectorTools::interpolate_boundary_values(
            *this->mapping,
            this->dof_handler,
            this->simulation_parameters.boundary_conditions.id[i_bc],
            dealii::Functions::ZeroFunction<dim>(dim + 1),
            this->zero_constraints,
            this->fe->component_mask(velocities));
        }
    }

  this->establish_solid_domain(false);

  this->zero_constraints.close();
}


template <int dim>
void
GLSNavierStokesSolver<dim>::setup_assemblers()
{
  this->assemblers.clear();

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

  // Buoyant force
  if (this->simulation_parameters.multiphysics.buoyancy_force)
    {
      this->assemblers.push_back(
        std::make_shared<BuoyancyAssembly<dim>>(this->simulation_control));
    }

  if (this->simulation_parameters.multiphysics.VOF)
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerBDF<dim>>(
              this->simulation_control));
        }

      // Surface tension force (STF)
      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          if (this->simulation_parameters.multiphysics.vof_parameters
                .surface_tension_force.enable_marangoni_effect)
            this->assemblers.push_back(
              std::make_shared<GLSNavierStokesVOFAssemblerMarangoni<dim>>(
                this->simulation_control,
                this->simulation_parameters.multiphysics.vof_parameters
                  .surface_tension_force));
          else
            {
              this->assemblers.push_back(
                std::make_shared<GLSNavierStokesVOFAssemblerSTF<dim>>(
                  this->simulation_control, this->simulation_parameters));
            }
        }

      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerNonNewtonianCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
      else
        {
          // Core assembler
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesVOFAssemblerCore<dim>>(
              this->simulation_control, this->simulation_parameters));
        }
    }
  else
    {
      // Time-stepping schemes
      if (is_bdf(this->simulation_control->get_assembly_method()))
        {
          if (!this->simulation_parameters.physical_properties_manager
                 .density_is_constant())
            {
              this->assemblers.push_back(
                std::make_shared<
                  GLSIsothermalCompressibleNavierStokesAssemblerBDF<dim>>(
                  this->simulation_control));
            }
          else
            {
              this->assemblers.push_back(
                std::make_shared<GLSNavierStokesAssemblerBDF<dim>>(
                  this->simulation_control));
            }
        }
      else if (is_sdirk(this->simulation_control->get_assembly_method()))
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSDIRK<dim>>(
              this->simulation_control));
        }

      // Velocity sources term
      if (this->simulation_parameters.velocity_sources.type ==
          Parameters::VelocitySource::VelocitySourceType::srf)
        {
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerSRF<dim>>(
              this->simulation_parameters.velocity_sources));
        }

      if (this->simulation_parameters.physical_properties_manager
            .is_non_newtonian())
        {
          // Core assembler with Non newtonian viscosity
          this->assemblers.push_back(
            std::make_shared<GLSNavierStokesAssemblerNonNewtonianCore<dim>>(
              this->simulation_control));
        }
      else
        {
          // Core default assembler
          if ((this->simulation_parameters.stabilization
                 .use_default_stabilization == true) ||
              this->simulation_parameters.stabilization.stabilization ==
                Parameters::Stabilization::NavierStokesStabilization::pspg_supg)
            {
              if (!this->simulation_parameters.physical_properties_manager
                     .density_is_constant())
                {
                  this->assemblers.push_back(
                    std::make_shared<
                      GLSIsothermalCompressibleNavierStokesAssemblerCore<dim>>(
                      this->simulation_control));
                }
              else
                {
                  this->assemblers.push_back(
                    std::make_shared<PSPGSUPGNavierStokesAssemblerCore<dim>>(
                      this->simulation_control));
                }
            }
          else if (this->simulation_parameters.stabilization.stabilization ==
                   Parameters::Stabilization::NavierStokesStabilization::gls)
            {
              if (!this->simulation_parameters.physical_properties_manager
                     .density_is_constant())
                {
                  this->assemblers.push_back(
                    std::make_shared<
                      GLSIsothermalCompressibleNavierStokesAssemblerCore<dim>>(
                      this->simulation_control));
                }
              else
                {
                  this->assemblers.push_back(
                    std::make_shared<GLSNavierStokesAssemblerCore<dim>>(
                      this->simulation_control));
                }
            }
          else
            throw std::runtime_error(
              "Using the GLS solver with a stabilization other than the pspg_supg or gls "
              "stabilization will lead to an unstable solver that is unable to converge");
        }
    }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_system_matrix()
{
  this->GLSNavierStokesSolver<
    dim>::assemble_system_matrix_without_preconditioner();

  // Assemble the preconditioner
  this->setup_preconditioner();
}



template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_system_matrix_without_preconditioner()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble matrix");

  this->system_matrix = 0;
  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          const DoFHandler<dim> *projected_phase_fraction_gradient_dof_handler =
            this->multiphysics
              ->get_projected_phase_fraction_gradient_dof_handler();
          const DoFHandler<dim> *curvature_dof_handler =
            this->multiphysics->get_curvature_dof_handler();

          scratch_data.enable_projected_phase_fraction_gradient(
            projected_phase_fraction_gradient_dof_handler->get_fe(),
            *this->cell_quadrature,
            *this->mapping);
          scratch_data.enable_curvature(curvature_dof_handler->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
        }
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }

  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &GLSNavierStokesSolver::assemble_local_system_matrix,
    &GLSNavierStokesSolver::copy_local_matrix_to_global_matrix,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));
  system_matrix.compress(VectorOperation::add);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->solution_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF),
        std::vector<TrilinosWrappers::MPI::Vector>());

      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          const DoFHandler<dim> *projected_phase_fraction_gradient_dof_handler =
            this->multiphysics
              ->get_projected_phase_fraction_gradient_dof_handler();

          typename DoFHandler<dim>::active_cell_iterator
            projected_phase_fraction_gradient_cell(
              &(*(this->triangulation)),
              cell->level(),
              cell->index(),
              projected_phase_fraction_gradient_dof_handler);
          scratch_data.reinit_projected_phase_fraction_gradient(
            projected_phase_fraction_gradient_cell,
            *this->multiphysics
               ->get_projected_phase_fraction_gradient_solution());



          const DoFHandler<dim> *curvature_dof_handler =
            this->multiphysics->get_curvature_dof_handler();
          typename DoFHandler<dim>::active_cell_iterator curvature_cell(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            curvature_dof_handler);
          scratch_data.reinit_curvature(
            curvature_cell, *this->multiphysics->get_curvature_solution());
        }
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      typename DoFHandler<dim>::active_cell_iterator temperature_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_ht);

      scratch_data.reinit_heat_transfer(temperature_cell,
                                        *this->multiphysics->get_solution(
                                          PhysicsID::heat_transfer));
    }

  scratch_data.calculate_physical_properties();

  copy_data.reset();

  if (this->simulation_parameters.physical_properties_manager
          .get_number_of_solids() < 1 ||
      cell->material_id() == 0)
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_matrix(scratch_data, copy_data);
        }
    }

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::copy_local_matrix_to_global_matrix(
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
GLSNavierStokesSolver<dim>::assemble_system_rhs()
{
  TimerOutput::Scope t(this->computing_timer, "Assemble RHS");

  this->system_rhs = 0;
  setup_assemblers();

  auto scratch_data = NavierStokesScratchData<dim>(
    this->simulation_parameters.physical_properties_manager,
    *this->fe,
    *this->cell_quadrature,
    *this->mapping,
    *this->face_quadrature);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      scratch_data.enable_vof(
        dof_handler_vof->get_fe(),
        *this->cell_quadrature,
        *this->mapping,
        this->simulation_parameters.multiphysics.vof_parameters.phase_filter);

      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          const DoFHandler<dim> *projected_phase_fraction_gradient_dof_handler =
            this->multiphysics
              ->get_projected_phase_fraction_gradient_dof_handler();
          const DoFHandler<dim> *curvature_dof_handler =
            this->multiphysics->get_curvature_dof_handler();

          scratch_data.enable_projected_phase_fraction_gradient(
            projected_phase_fraction_gradient_dof_handler->get_fe(),
            *this->cell_quadrature,
            *this->mapping);
          scratch_data.enable_curvature(curvature_dof_handler->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
        }
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);
      scratch_data.enable_heat_transfer(dof_handler_ht->get_fe(),
                                        *this->cell_quadrature,
                                        *this->mapping);
    }

  WorkStream::run(
    this->dof_handler.begin_active(),
    this->dof_handler.end(),
    *this,
    &GLSNavierStokesSolver::assemble_local_system_rhs,
    &GLSNavierStokesSolver::copy_local_rhs_to_global_rhs,
    scratch_data,
    StabilizedMethodsTensorCopyData<dim>(this->fe->n_dofs_per_cell(),
                                         this->cell_quadrature->size()));

  this->system_rhs.compress(VectorOperation::add);

  if (this->simulation_control->is_first_assembly())
    this->simulation_control->provide_residual(this->system_rhs.l2_norm());
}


template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  NavierStokesScratchData<dim> &                        scratch_data,
  StabilizedMethodsTensorCopyData<dim> &                copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(
    cell,
    this->evaluation_point,
    this->previous_solutions,
    this->solution_stages,
    this->forcing_function,
    this->flow_control.get_beta(),
    this->simulation_parameters.stabilization.pressure_scaling_factor);

  if (this->simulation_parameters.multiphysics.VOF)
    {
      const DoFHandler<dim> *dof_handler_vof =
        this->multiphysics->get_dof_handler(PhysicsID::VOF);
      typename DoFHandler<dim>::active_cell_iterator phase_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_vof);

      scratch_data.reinit_vof(
        phase_cell,
        *this->multiphysics->get_solution(PhysicsID::VOF),
        *this->multiphysics->get_filtered_solution(PhysicsID::VOF),
        *this->multiphysics->get_previous_solutions(PhysicsID::VOF),
        std::vector<TrilinosWrappers::MPI::Vector>());

      if (this->simulation_parameters.multiphysics.vof_parameters
            .surface_tension_force.enable)
        {
          const DoFHandler<dim> *projected_phase_fraction_gradient_dof_handler =
            this->multiphysics
              ->get_projected_phase_fraction_gradient_dof_handler();
          typename DoFHandler<dim>::active_cell_iterator
            projected_phase_fraction_gradient_cell(
              &(*(this->triangulation)),
              cell->level(),
              cell->index(),
              projected_phase_fraction_gradient_dof_handler);
          scratch_data.reinit_projected_phase_fraction_gradient(
            projected_phase_fraction_gradient_cell,
            *this->multiphysics
               ->get_projected_phase_fraction_gradient_solution());

          const DoFHandler<dim> *curvature_dof_handler =
            this->multiphysics->get_curvature_dof_handler();
          typename DoFHandler<dim>::active_cell_iterator curvature_cell(
            &(*(this->triangulation)),
            cell->level(),
            cell->index(),
            curvature_dof_handler);
          scratch_data.reinit_curvature(
            curvature_cell, *this->multiphysics->get_curvature_solution());
        }
    }

  if (this->simulation_parameters.multiphysics.heat_transfer)
    {
      const DoFHandler<dim> *dof_handler_ht =
        this->multiphysics->get_dof_handler(PhysicsID::heat_transfer);

      typename DoFHandler<dim>::active_cell_iterator temperature_cell(
        &(*(this->triangulation)),
        cell->level(),
        cell->index(),
        dof_handler_ht);

      scratch_data.reinit_heat_transfer(temperature_cell,
                                        *this->multiphysics->get_solution(
                                          PhysicsID::heat_transfer));
    }

  scratch_data.calculate_physical_properties();

  copy_data.reset();

  if (this->simulation_parameters.physical_properties_manager
          .get_number_of_solids() < 1 ||
      cell->material_id() == 0)
    {
      for (auto &assembler : this->assemblers)
        {
          assembler->assemble_rhs(scratch_data, copy_data);
        }
    }


  cell->get_dof_indices(copy_data.local_dof_indices);
}



template <int dim>
void
GLSNavierStokesSolver<dim>::copy_local_rhs_to_global_rhs(
  const StabilizedMethodsTensorCopyData<dim> &copy_data)
{
  if (!copy_data.cell_is_local)
    return;

  const AffineConstraints<double> &constraints_used = this->zero_constraints;
  constraints_used.distribute_local_to_global(copy_data.local_rhs,
                                              copy_data.local_dof_indices,
                                              this->system_rhs);
}

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
GLSNavierStokesSolver<dim>::set_initial_condition_fd(
  Parameters::InitialConditionType initial_condition_type,
  bool                             restart)
{
  if (restart)
    {
      this->pcout << "************************" << std::endl;
      this->pcout << "---> Simulation Restart " << std::endl;
      this->pcout << "************************" << std::endl;
      this->read_checkpoint();
    }
  else if (initial_condition_type ==
           Parameters::InitialConditionType::L2projection)
    {
      assemble_L2_projection();
      setup_preconditioner();

      if (this->simulation_parameters.linear_solver.solver ==
          Parameters::LinearSolver::SolverType::amg)
        solve_system_AMG(true, 1e-15, 1e-15);
      else
        solve_system_GMRES(true, 1e-15, 1e-15);

      this->present_solution = this->newton_update;
      this->finish_time_step();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::nodal)
    {
      this->set_nodal_values();
      this->finish_time_step();
    }
  else if (initial_condition_type == Parameters::InitialConditionType::viscous)
    {
      this->set_nodal_values();
      std::shared_ptr<RheologicalModel> original_viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Temporarily set the rheology to be newtonian with predefined viscosity
      std::shared_ptr<Newtonian> temporary_rheology =
        std::make_shared<Newtonian>(
          this->simulation_parameters.initial_condition->viscosity);

      this->simulation_parameters.physical_properties_manager.set_rheology(
        temporary_rheology);


      this->simulation_control->set_assembly_method(
        Parameters::SimulationControl::TimeSteppingMethod::steady);
      PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
        false);
      this->finish_time_step();

      this->simulation_parameters.physical_properties_manager.set_rheology(
        original_viscosity_model);
    }
  else if (initial_condition_type == Parameters::InitialConditionType::ramp)
    {
      this->pcout << "*********************************" << std::endl;
      this->pcout << " Initial condition using ramp " << std::endl;
      this->pcout << "*********************************" << std::endl;

      this->set_nodal_values();

      // Create a pointer to the current viscosity model
      std::shared_ptr<RheologicalModel> viscosity_model =
        this->simulation_parameters.physical_properties_manager.get_rheology();

      // Gather n and viscosity final paramerters
      const double n_end         = viscosity_model->get_n();
      const double viscosity_end = viscosity_model->get_viscosity();

      // n ramp parameters
      double n =
        this->simulation_parameters.initial_condition->ramp.ramp_n.n_init;
      const int n_iter_n =
        this->simulation_parameters.initial_condition->ramp.ramp_n.n_iter;
      const double alpha_n =
        this->simulation_parameters.initial_condition->ramp.ramp_n.alpha;

      // viscosity ramp parameters
      const int n_iter_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .n_iter;
      double       viscosity = n_iter_viscosity > 0 ?
                                 this->simulation_parameters.initial_condition->ramp
                             .ramp_viscosity.viscosity_init :
                                 viscosity_end;
      const double alpha_viscosity =
        this->simulation_parameters.initial_condition->ramp.ramp_viscosity
          .alpha;

      viscosity_model->set_viscosity(viscosity);

      // Ramp on n
      for (int i = 0; i < n_iter_n; ++i)
        {
          this->pcout << std::setprecision(4)
                      << "********* Solution for n = " + std::to_string(n) +
                           " and viscosity = " + std::to_string(viscosity) +
                           " *********"
                      << std::endl;

          viscosity_model->set_n(n);

          this->simulation_control->set_assembly_method(
            Parameters::SimulationControl::TimeSteppingMethod::steady);
          PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
            false);
          this->finish_time_step();

          n += alpha_n * (n_end - n);
        }

      // Reset n to simulation parameters
      viscosity_model->set_n(n_end);

      // Ramp on viscosity
      for (int i = 0; i < n_iter_viscosity; ++i)
        {
          this->pcout << std::setprecision(4)
                      << "********* Solution for n = " + std::to_string(n_end) +
                           " and viscosity = " + std::to_string(viscosity) +
                           " *********"
                      << std::endl;

          viscosity_model->set_viscosity(viscosity);

          this->simulation_control->set_assembly_method(
            Parameters::SimulationControl::TimeSteppingMethod::steady);
          PhysicsSolver<TrilinosWrappers::MPI::Vector>::solve_non_linear_system(
            false);
          this->finish_time_step();

          viscosity += alpha_viscosity * (viscosity_end - viscosity);
        }

      // Reset viscosity to simulation parameters
      viscosity_model->set_viscosity(viscosity_end);
    }
  else
    {
      throw std::runtime_error("GLSNS - Initial condition could not be set");
    }
}

template <int dim>
void
GLSNavierStokesSolver<dim>::assemble_L2_projection()
{
  system_matrix    = 0;
  this->system_rhs = 0;
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

  for (const auto &cell : this->dof_handler.active_cell_iterators())
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
            this->system_rhs);
        }
    }
  system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_linear_system(const bool initial_step,
                                                const bool /* renewed_matrix */)
{
  const double absolute_residual =
    this->simulation_parameters.linear_solver.minimum_residual;
  const double relative_residual =
    this->simulation_parameters.linear_solver.relative_residual;

  if (this->simulation_parameters.linear_solver.solver ==
      Parameters::LinearSolver::SolverType::gmres)
    solve_system_GMRES(initial_step, absolute_residual, relative_residual);
  else if (this->simulation_parameters.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::bicgstab)
    solve_system_BiCGStab(initial_step, absolute_residual, relative_residual);
  else if (this->simulation_parameters.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::amg)
    solve_system_AMG(initial_step, absolute_residual, relative_residual);
  else if (this->simulation_parameters.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::direct)
    solve_system_direct(initial_step, absolute_residual, relative_residual);
  else
    throw(std::runtime_error("This solver is not allowed"));
  this->rescale_pressure_dofs_in_newton_update();
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_preconditioner()
{
  if (this->simulation_parameters.linear_solver.solver ==
        Parameters::LinearSolver::SolverType::gmres ||
      this->simulation_parameters.linear_solver.solver ==
        Parameters::LinearSolver::SolverType::bicgstab)
    setup_ILU();
  else if (this->simulation_parameters.linear_solver.solver ==
           Parameters::LinearSolver::SolverType::amg)
    setup_AMG();
}


template <int dim>
void
GLSNavierStokesSolver<dim>::setup_ILU()
{
  TimerOutput::Scope t(this->computing_timer, "setup_ILU");

  const double ilu_atol =
    this->simulation_parameters.linear_solver.ilu_precond_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.ilu_precond_rtol;
  TrilinosWrappers::PreconditionILU::AdditionalData preconditionerOptions(
    current_preconditioner_fill_level, ilu_atol, ilu_rtol, 0);

  ilu_preconditioner = std::make_shared<TrilinosWrappers::PreconditionILU>();

  ilu_preconditioner->initialize(system_matrix, preconditionerOptions);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::setup_AMG()
{
  TimerOutput::Scope t(this->computing_timer, "setup_AMG");

  // Constant modes for velocity
  std::vector<std::vector<bool>> constant_modes;

  // Constant modes include pressure since everything is in the same matrix
  ComponentMask components(dim + 1, true);
  DoFTools::extract_constant_modes(this->dof_handler,
                                   components,
                                   constant_modes);

  TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
  amg_data.constant_modes = constant_modes;

  const bool elliptic              = false;
  bool       higher_order_elements = false;
  if (this->velocity_fem_degree > 1)
    higher_order_elements = true;
  const unsigned int n_cycles =
    this->simulation_parameters.linear_solver.amg_n_cycles;
  const bool   w_cycle = this->simulation_parameters.linear_solver.amg_w_cycles;
  const double aggregation_threshold =
    this->simulation_parameters.linear_solver.amg_aggregation_threshold;
  const unsigned int smoother_sweeps =
    this->simulation_parameters.linear_solver.amg_smoother_sweeps;
  const unsigned int smoother_overlap =
    this->simulation_parameters.linear_solver.amg_smoother_overlap;
  const bool                                        output_details = false;
  const char *                                      smoother_type  = "ILU";
  const char *                                      coarse_type    = "ILU";
  TrilinosWrappers::PreconditionAMG::AdditionalData preconditionerOptions(
    elliptic,
    higher_order_elements,
    n_cycles,
    w_cycle,
    aggregation_threshold,
    constant_modes,
    smoother_sweeps,
    smoother_overlap,
    output_details,
    smoother_type,
    coarse_type);

  Teuchos::ParameterList              parameter_ml;
  std::unique_ptr<Epetra_MultiVector> distributed_constant_modes;
  preconditionerOptions.set_parameters(parameter_ml,
                                       distributed_constant_modes,
                                       system_matrix);
  const double ilu_fill = current_preconditioner_fill_level;
  const double ilu_atol =
    this->simulation_parameters.linear_solver.amg_precond_ilu_atol;
  const double ilu_rtol =
    this->simulation_parameters.linear_solver.amg_precond_ilu_rtol;
  parameter_ml.set("smoother: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("smoother: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("smoother: ifpack relative threshold", ilu_rtol);

  parameter_ml.set("coarse: ifpack level-of-fill", ilu_fill);
  parameter_ml.set("coarse: ifpack absolute threshold", ilu_atol);
  parameter_ml.set("coarse: ifpack relative threshold", ilu_rtol);
  amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>();
  amg_preconditioner->initialize(system_matrix, parameter_ml);
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_GMRES(const bool   initial_step,
                                               const double absolute_residual,
                                               const double relative_residual)
{
  // Try multiple fill of the ILU preconditioner. Start from the initial fill
  // given in the parameter file. If for any reason the linear solver would have
  // crash it will restart with a fill level increased by 1. This restart
  // process will happen up to a maximum of 20 times, after which it will let
  // the solver crash. if a change happened on the fill level it will go back to
  // its original value at the end of the restart process.

  const unsigned int max_iter = 3;
  unsigned int       iter     = 0;
  bool               success  = false;


  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);
  bool extra_verbose = false;
  if (this->simulation_parameters.linear_solver.verbosity ==
      Parameters::Verbosity::extra_verbose)
    extra_verbose = true;

  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    extra_verbose,
    this->simulation_parameters.linear_solver.max_krylov_vectors);
  while (success == false and iter < max_iter)
    {
      try
        {
          if (!ilu_preconditioner)
            setup_preconditioner();

          TrilinosWrappers::SolverGMRES solver(solver_control,
                                               solver_parameters);

          {
            TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

            solver.solve(system_matrix,
                         completely_distributed_solution,
                         system_rhs,
                         *ilu_preconditioner);

            if (this->simulation_parameters.linear_solver.verbosity !=
                Parameters::Verbosity::quiet)
              {
                this->pcout
                  << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
              }
          }
          constraints_used.distribute(completely_distributed_solution);
          auto &newton_update = this->newton_update;
          newton_update       = completely_distributed_solution;
          success             = true;
        }
      catch (std::exception &e)
        {
          current_preconditioner_fill_level += 1;
          this->pcout
            << " GMRES solver failed! Trying with a higher preconditioner fill level. New fill = "
            << current_preconditioner_fill_level << std::endl;
          setup_preconditioner();

          if (iter == max_iter - 1 && !this->simulation_parameters.linear_solver
                                         .force_linear_solver_continuation)
            throw e;
        }
      iter += 1;
    }
  current_preconditioner_fill_level = initial_preconditioner_fill_level;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_BiCGStab(
  const bool   initial_step,
  const double absolute_residual,
  const double relative_residual)
{
  TimerOutput::Scope t(this->computing_timer, "solve");
  const unsigned int max_iter = 3;
  unsigned int       iter     = 0;
  bool               success  = false;


  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);
  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  bool extra_verbose = false;
  if (this->simulation_parameters.linear_solver.verbosity ==
      Parameters::Verbosity::extra_verbose)
    extra_verbose = true;
  TrilinosWrappers::SolverBicgstab::AdditionalData solver_parameters(
    extra_verbose);

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);
  TrilinosWrappers::SolverBicgstab solver(solver_control, solver_parameters);
  while (success == false and iter < max_iter)
    {
      try
        {
          if (!ilu_preconditioner)
            setup_preconditioner();

          {
            TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

            solver.solve(system_matrix,
                         completely_distributed_solution,
                         system_rhs,
                         *ilu_preconditioner);

            if (this->simulation_parameters.linear_solver.verbosity !=
                Parameters::Verbosity::quiet)
              {
                this->pcout
                  << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
              }
            constraints_used.distribute(completely_distributed_solution);
            this->newton_update = completely_distributed_solution;
          }
          success = true;
        }
      catch (std::exception &e)
        {
          current_preconditioner_fill_level += 1;
          this->pcout
            << " BiCGStab solver failed! Trying with a higher preconditioner fill level. New fill = "
            << current_preconditioner_fill_level << std::endl;
          setup_preconditioner();

          if (iter == max_iter - 1 && !this->simulation_parameters.linear_solver
                                         .force_linear_solver_continuation)
            throw e;
        }
      iter += 1;
    }
  current_preconditioner_fill_level = initial_preconditioner_fill_level;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_AMG(const bool   initial_step,
                                             const double absolute_residual,
                                             const double relative_residual)
{
  const unsigned int max_iter = 3;
  unsigned int       iter     = 0;
  bool               success  = false;


  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;

  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);
  if (this->simulation_parameters.linear_solver.verbosity !=
      Parameters::Verbosity::quiet)
    {
      this->pcout << "  -Tolerance of iterative solver is : "
                  << linear_solver_tolerance << std::endl;
    }
  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);
  bool extra_verbose = false;
  if (this->simulation_parameters.linear_solver.verbosity ==
      Parameters::Verbosity::extra_verbose)
    extra_verbose = true;
  TrilinosWrappers::SolverGMRES::AdditionalData solver_parameters(
    extra_verbose,
    this->simulation_parameters.linear_solver.max_krylov_vectors);

  TrilinosWrappers::SolverGMRES solver(solver_control, solver_parameters);
  while (success == false and iter < max_iter)
    {
      try
        {
          if (!amg_preconditioner)
            setup_preconditioner();

          {
            TimerOutput::Scope t(this->computing_timer, "solve_linear_system");

            solver.solve(system_matrix,
                         completely_distributed_solution,
                         system_rhs,
                         *amg_preconditioner);

            if (this->simulation_parameters.linear_solver.verbosity !=
                Parameters::Verbosity::quiet)
              {
                this->pcout
                  << "  -Iterative solver took : " << solver_control.last_step()
                  << " steps " << std::endl;
              }

            constraints_used.distribute(completely_distributed_solution);

            this->newton_update = completely_distributed_solution;
            success             = true;
          }
        }
      catch (std::exception &e)
        {
          current_preconditioner_fill_level += 1;
          this->pcout
            << " AMG solver failed! Trying with a higher preconditioner fill level. New fill = "
            << current_preconditioner_fill_level << std::endl;
          setup_preconditioner();

          if (iter == max_iter - 1 && !this->simulation_parameters.linear_solver
                                         .force_linear_solver_continuation)
            throw e;
        }
      iter += 1;
    }
  current_preconditioner_fill_level = initial_preconditioner_fill_level;
}


template <int dim>
void
GLSNavierStokesSolver<dim>::solve_system_direct(const bool   initial_step,
                                                const double absolute_residual,
                                                const double relative_residual)
{
  auto &system_rhs          = this->system_rhs;
  auto &nonzero_constraints = this->nonzero_constraints;

  const AffineConstraints<double> &constraints_used =
    initial_step ? nonzero_constraints : this->zero_constraints;
  const double linear_solver_tolerance =
    std::max(relative_residual * system_rhs.l2_norm(), absolute_residual);

  TrilinosWrappers::MPI::Vector completely_distributed_solution(
    this->locally_owned_dofs, this->mpi_communicator);

  SolverControl solver_control(
    this->simulation_parameters.linear_solver.max_iterations,
    linear_solver_tolerance,
    true,
    true);
  TrilinosWrappers::SolverDirect solver(solver_control);

  solver.initialize(system_matrix);
  solver.solve(completely_distributed_solution, system_rhs);
  constraints_used.distribute(completely_distributed_solution);
  auto &newton_update = this->newton_update;
  newton_update       = completely_distributed_solution;
}

template <int dim>
void
GLSNavierStokesSolver<dim>::solve()
{
  // This is enforced to 1 right now because it does not provide
  // better speed-up than using MPI. This could be eventually changed...
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();
  this->box_refine_mesh();
  this->set_initial_condition(
    this->simulation_parameters.initial_condition->type,
    this->simulation_parameters.restart_parameters.restart);
  this->update_multiphysics_time_average_solution();

  while (this->simulation_control->integrate())
    {
      this->forcing_function->set_time(
        this->simulation_control->get_current_time());

      if ((this->simulation_control->get_step_number() %
               this->simulation_parameters.mesh_adaptation.frequency !=
             0 ||
           this->simulation_parameters.mesh_adaptation.type ==
             Parameters::MeshAdaptation::Type::none ||
           this->simulation_control->is_at_start()) &&
          this->simulation_parameters.boundary_conditions.time_dependent)
        {
          update_boundary_conditions();
        }

      this->simulation_control->print_progression(this->pcout);
      this->dynamic_flow_control();

      if (this->simulation_control->is_at_start())
        {
          this->iterate();
        }
      else
        {
          NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>::
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
template class GLSNavierStokesSolver<2>;
template class GLSNavierStokesSolver<3>;
