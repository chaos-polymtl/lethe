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
  : NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>(p_nsparam)
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
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::setup_dofs_fd()
{
  // TODO
}

template <int dim>
void
MFGLSNavierStokesSolver<dim>::update_boundary_conditions()
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
MFGLSNavierStokesSolver<dim>::setup_operator()
{}

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
