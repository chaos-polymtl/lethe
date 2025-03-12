
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

#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>

template <int dim>
FluidDynamicsVANSMatrixFree<dim>::FluidDynamicsVANSMatrixFree(
  CFDDEMSimulationParameters<dim> &param)
  : FluidDynamicsMatrixFree<dim>(param.cfd_parameters)
  , cfd_dem_simulation_parameters(param)
  , particle_mapping(1)
  , particle_handler(
      *this->triangulation,
      particle_mapping,
      DEM::get_number_properties<DEM::CFDDEMProperties::PropertiesIndex>())
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
          mf_operator->evaluate_void_fraction(void_fraction_manager);
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

// Pre-compile the 2D and 3D solver to ensure that the
// library is valid before we actually compile the solver
template class FluidDynamicsVANSMatrixFree<2>;
template class FluidDynamicsVANSMatrixFree<3>;
