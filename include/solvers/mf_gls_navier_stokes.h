/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2023 by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------*/

#ifndef lethe_mf_gls_navier_stokes_h
#define lethe_mf_gls_navier_stokes_h

#include <core/exceptions.h>

#include <solvers/navier_stokes_base.h>

using namespace dealii;

/**
 * @brief A solver class for the Navier-Stokes equation using GLS stabilization and the matrix free approach
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 */

template <int dim>
class MFGLSNavierStokesSolver
  : public NavierStokesBase<dim,
                            LinearAlgebra::distributed::Vector<double>,
                            IndexSet>
{
public:
  MFGLSNavierStokesSolver(SimulationParameters<dim> &nsparam);
  ~MFGLSNavierStokesSolver();

  /**
   * @brief solve Solves the Navier-Stokes problem
   *
   * This functions solves the problem defined by the Navier-Stokes paramter
   * by iterating through time or through the mesh refinements.
   */
  virtual void
  solve();

protected:
  /**
   * @brief setup_dofs_fd Setup the degree of freedom, system matrix and solution vectors
   * Navier-Stokes problem
   */
  virtual void
  setup_dofs_fd() override;

  /**
   * @brief update non zero constraint if the boundary is time dependent
   */
  void
  update_boundary_conditions();

  /**
   * @brief Sets the initial condition for the solver
   *
   * If the simulation is restarted from a checkpoint, the initial solution
   * setting is bypassed and the checkpoint is instead read.
   *
   * @param initial_condition_type The type of initial condition to be set
   *
   * @param restart A boolean that indicates if the simulation is being restarted.
   * if set to true, the initial conditions are never set, but are instead
   * overriden by the read_checkpoint functionnality.
   *
   **/
  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) override;

protected:
  /**
   *  @brief Assembles the matrix associated with the solver
   */
  virtual void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the rhs associated with the solver
   */
  virtual void
  assemble_system_rhs() override;

  /**
   * @brief  Updates the average velocity field solution in the multiphyscics interface
   */
  virtual void
  update_multiphysics_time_average_solution() override;

  /**
   * @brief  Set-up the appropriate preconditioner
   */
  void
  setup_preconditioner();

  /**
   * @brief  defined the non zero constraints used to solved the problem.
   */
  void
  define_non_zero_constraints();

  /**
   * @brief defined the zero_constraints used to solved the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief sets up the operator
   */
  virtual void
  setup_operator();

  /**
   * @brief Call for the assembly of the linear system of equation
   *
   * @param initial_step Indicates if this is the first solution of the linear system.
   * If this is the case, the non_zero version of the constraints are used for
   * the Dirichlet boundary conditions
   *
   * @param renewed_matrix Indicates if the matrix has been reassembled, and thus
   * the preconditioner needs to be reassmbled.
   *
   * //TODO the renewed_matrix parameters needs to be deprecated
   *
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

private:
  /**
   * @brief Assembles an L2_projection matrix for the velocity and the pressure.
   * This L2 projection matrix can be used to set the initial condition for
   * the Navier-Stokes problem
   */
  void
  assemble_L2_projection();

  /**
   * @brief GMRES solver with GMG preconditioning
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual);

  /**
   * @brief  Set-up GMG preconditioner
   */
  void
  setup_GMG();

  /**
   * Members
   */
protected:
  // TODO create operator: NSOperator<dim, double> system_matrix;
};

#endif
