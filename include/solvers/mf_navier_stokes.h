/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
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

#ifndef lethe_mf_navier_stokes_h
#define lethe_mf_navier_stokes_h

#include <core/exceptions.h>

#include <solvers/mf_navier_stokes_operators.h>
#include <solvers/navier_stokes_base.h>

using namespace dealii;

/**
 * @brief A solver class for the Navier-Stokes equation that uses the matrix
 * free approach.
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved.
 */
template <int dim>
class MFNavierStokesSolver
  : public NavierStokesBase<dim,
                            LinearAlgebra::distributed::Vector<double>,
                            IndexSet>
{
  using VectorType = LinearAlgebra::distributed::Vector<double>;

public:
  /**
   * @brief Construct a new MFNavierStokesSolver object.
   *
   * @param nsparam Relevant parameters for the solver.
   */
  MFNavierStokesSolver(SimulationParameters<dim> &nsparam);

  /**
   * @brief Destroy the MFNavierStokesSolver object.
   *
   */
  ~MFNavierStokesSolver();

  /**
   * @brief Solves the Navier-Stokes problem.
   *
   * This functions solves the problem defined by the Navier-Stokes parameter
   * file by iterating through time or through the mesh refinements.
   */
  virtual void
  solve();

protected:
  /**
   * @brief Setup the DoFs, system operator and solution vectors.
   */
  virtual void
  setup_dofs_fd() override;

  /**
   * @brief Update non zero constraints if the boundary is time dependent.
   */
  void
  update_boundary_conditions();

  /**
   * @brief Set the initial conditions for the solver.
   *
   * If the simulation is restarted from a checkpoint, the initial solution
   * setting is bypassed and the checkpoint is instead read.
   *
   * @param initial_condition_type The type of initial condition to be set.
   *
   * @param restart A boolean that indicates if the simulation is being restarted.
   * If set to true, the initial conditions are never set, but are instead
   * overriden by the read_checkpoint functionality.
   *
   **/
  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) override;

  /**
   * @brief Assemble the matrix associated with the solver. This function is only
   * required for compilation and it is not used for matrix free solvers.
   */
  virtual void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the system rhs associated with the solver.
   */
  virtual void
  assemble_system_rhs() override;

  /**
   * @brief  Update the average velocity field solution in the multiphyscics interface.
   */
  virtual void
  update_multiphysics_time_average_solution() override;

  /**
   * @brief  Set up the appropriate preconditioner.
   */
  void
  setup_preconditioner();

  /**
   * @brief Define the non-zero constraints used to solve the problem.
   */
  void
  define_non_zero_constraints();

  /**
   * @brief Define the zero constraints used to solve the problem.
   */
  void
  define_zero_constraints();

  /**
   * @brief Call for the assembly of the linear system of equations.
   *
   * @param initial_step Indicates if this is the first solution of the linear system.
   * If this is the case, the non_zero version of the constraints are used for
   * the Dirichlet boundary conditions.
   *
   * @param renewed_matrix Indicates if the matrix has been reassembled, and thus
   * the preconditioner needs to be reassmbled.
   *
   * TODO: the renewed_matrix parameter needs to be deprecated.
   *
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

private:
  /**
   * @brief Assembles an L2_projection matrix for the velocity and the pressure.
   * This L2 projection matrix can be used to set the initial condition for
   * the Navier-Stokes problem.
   */
  void
  assemble_L2_projection();

  /**
   * @brief GMRES solver with preconditioning.
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual);

protected:
  // Matrix-free operator
  NavierStokesOperatorBase<dim, double> *system_operator;
};

#endif
