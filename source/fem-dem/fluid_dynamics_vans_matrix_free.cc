
// SPDX-FileCopyrightText: Copyright (c) 2023-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/grids.h>
#include <core/manifolds.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <fem-dem/fluid_dynamics_vans_matrix_free.h>
#include <fem-dem/fluid_dynamics_vans_matrix_free_operators.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>


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
  const VoidFractionBase<dim>              &void_fraction_manager)
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

      this->void_fraction_dof_handlers.resize(min_level, max_level);

      for (unsigned int l = min_level; l <= max_level; l++)
        {
          this->void_fraction_dof_handlers[l].reinit(
            this->dof_handlers[l].get_triangulation());
          this->void_fraction_dof_handlers[l].distribute_dofs(
            void_fraction_manager.dof_handler.get_fe());
        }

      this->transfers_void_fraction.resize(min_level, max_level);

      for (unsigned int l = min_level; l < max_level; l++)
        {
          this->transfers_void_fraction[l + 1].reinit(
            this->void_fraction_dof_handlers[l + 1],
            this->void_fraction_dof_handlers[l],
            {},
            {});
        }

      this->mg_transfer_gc_void_fraction =
        std::make_shared<MFNavierStokesVANSPreconditionGMG::GCTransferType>(
          this->transfers_void_fraction);

#if DEAL_II_VERSION_GTE(9, 7, 0)
      this->mg_transfer_gc_void_fraction->build(
        void_fraction_manager.dof_handler, [&](const auto l, auto &vec) {
          vec.reinit(
            this->void_fraction_dof_handlers[l].locally_owned_dofs(),
            DoFTools::extract_locally_active_dofs(
              this->void_fraction_dof_handlers[l]),
            this->void_fraction_dof_handlers[l].get_mpi_communicator());
        });
#endif

      MGLevelObject<MFNavierStokesVANSPreconditionGMG::MGVectorType>
        mg_void_fraction_solution(this->minlevel, this->maxlevel);


      // A deal.II vector is required here, so we take the deal.II vector
      // solution from the void fraction manager instead of the trilinos vector
      // one.
      this->mg_transfer_gc_void_fraction->interpolate_to_mg(
        void_fraction_manager.dof_handler,
        mg_void_fraction_solution,
        void_fraction_manager.void_fraction_solution);

      for (unsigned int l = min_level; l <= max_level; l++)
        {
          mg_void_fraction_solution[l].update_ghost_values();

          if (auto mf_operator = dynamic_cast<VANSOperator<dim, double> *>(
                &(*this->mg_operators[l])))
            mf_operator->compute_void_fraction(
              mg_void_fraction_solution[l],
              this->void_fraction_dof_handlers[l]);
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
  , void_fraction_manager(
      &(*this->triangulation),
      param.void_fraction,
      this->cfd_dem_simulation_parameters.cfd_parameters.linear_solver.at(
        PhysicsID::fluid_dynamics),
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

  void_fraction_manager.setup_dofs();
  void_fraction_manager.setup_constraints(
    this->cfd_dem_simulation_parameters.cfd_parameters.boundary_conditions);
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::finish_time_step_fd()
{
  // Void fraction percolation must be done before the time step is finished to
  // ensure that the checkpointed information is correct
  void_fraction_manager.percolate_void_fraction();

  FluidDynamicsMatrixFree<dim>::finish_time_step();
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::output_field_hook(DataOut<dim> &data_out)
{
  data_out.add_data_vector(void_fraction_manager.dof_handler,
                           void_fraction_manager.void_fraction_locally_relevant,
                           "void_fraction");
}

template <int dim>
void
FluidDynamicsVANSMatrixFree<dim>::solve()
{
  this->computing_timer.enter_subsection("Read mesh and manifolds");

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->computing_timer.leave_subsection("Read mesh and manifolds");

  this->setup_dofs();

  void_fraction_manager.calculate_void_fraction(
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

      if (is_bdf(this->simulation_control->get_assembly_method()))
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

      // Calculate the void fraction and evaluate it within the matrix-free
      // operator
      {
        TimerOutput::Scope t(this->computing_timer, "Calculate void fraction");
        void_fraction_manager.calculate_void_fraction(
          this->simulation_control->get_current_time());

        // The base matrix-free operator is not aware of the void fraction. We
        // must do a cast here to ensure that the operator is of the right type
        if (auto mf_operator = dynamic_cast<VANSOperator<dim, double> *>(
              this->system_operator.get()))
          mf_operator->compute_void_fraction(
            void_fraction_manager.void_fraction_solution,
            void_fraction_manager.dof_handler);
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
      this->dof_handler,
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
                 this->present_solution,
                 this->time_derivative_previous_solutions,
                 this->void_fraction_manager);
}

// Pre-compile the 2D and 3D solver to ensure that the
// library is valid before we actually compile the solver
template class FluidDynamicsVANSMatrixFree<2>;
template class FluidDynamicsVANSMatrixFree<3>;
