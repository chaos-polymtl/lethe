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
 * ---------------------------------------------------------------------*/

#include "solvers/mf_gls_navier_stokes.h"

#include <core/bdf.h>
#include <core/grids.h>
#include <core/manifolds.h>
#include <core/multiphysics.h>
#include <core/sdirk.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>

#include <solvers/mf_gls_navier_stokes.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools.h>

// Constructor for class MFGLSNavierStokesSolver
template <int dim>
MFGLSNavierStokesSolver<dim>::MFGLSNavierStokesSolver(
  SimulationParameters<dim> &p_nsparam)
  : NavierStokesBase<dim, LinearAlgebra::distributed::Vector<double>, IndexSet>(
      p_nsparam)
{
  // TODO
}

template <int dim>
MFGLSNavierStokesSolver<dim>::~MFGLSNavierStokesSolver()
{
  this->dof_handler.clear();
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::solve()
{
  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();

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

      // if (this->simulation_control->is_at_start())
      //   {
      //     this->iterate();
      //   }
      // else
      //   {
      NavierStokesBase<dim,
                       LinearAlgebra::distributed::Vector<double>,
                       IndexSet>::refine_mesh();
      //     this->iterate();
      //   }

      this->postprocess(false);
      this->finish_time_step();
    }

  this->finish_simulation();
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::setup_dofs_fd()
{
  // Fill the dof handler and initialize vectors
  this->dof_handler.distribute_dofs(*this->fe);
  DoFRenumbering::Cuthill_McKee(this->dof_handler);

  this->locally_owned_dofs = this->dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs(this->dof_handler,
                                          this->locally_relevant_dofs);

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

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
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
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::update_boundary_conditions()
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

/**
 * Set the initial condition using a L2 or a viscous solver
 **/
template <int dim>
void
MFGLSNavierStokesSolver<dim>::set_initial_condition_fd(
  Parameters::InitialConditionType /* initial_condition_type */,
  bool /* restart */)
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::assemble_system_matrix()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::assemble_system_rhs()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::update_multiphysics_time_average_solution()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::setup_preconditioner()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::define_non_zero_constraints()
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
MFGLSNavierStokesSolver<dim>::define_zero_constraints()
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
MFGLSNavierStokesSolver<dim>::setup_operator()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::solve_linear_system(
  const bool /* initial_step */,
  const bool /* renewed_matrix */)
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::assemble_L2_projection()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::solve_system_GMRES(
  const bool /* initial_step */,
  const double /* absolute_residual */,
  const double /* relative_residual */)
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::setup_GMG()
{
  // TODO
}

// Pre-compile the 2D and 3D Navier-Stokes solver to ensure that the library is
// valid before we actually compile the solver This greatly helps with debugging
template class MFGLSNavierStokesSolver<2>;
template class MFGLSNavierStokesSolver<3>;
