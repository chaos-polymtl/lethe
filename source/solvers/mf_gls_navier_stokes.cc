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
  MultithreadInfo::set_thread_limit(1);

  read_mesh_and_manifolds(
    *this->triangulation,
    this->simulation_parameters.mesh,
    this->simulation_parameters.manifolds_parameters,
    this->simulation_parameters.restart_parameters.restart,
    this->simulation_parameters.boundary_conditions);

  this->setup_dofs();

  while (this->simulation_control->integrate())
    {
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

  this->newton_update.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->system_rhs.reinit(this->locally_owned_dofs, this->mpi_communicator);
  this->local_evaluation_point.reinit(this->locally_owned_dofs,
                                      this->mpi_communicator);
  auto &nonzero_constraints = this->get_nonzero_constraints();

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
  // TODO
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
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::define_zero_constraints()
{
  // TODO
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
