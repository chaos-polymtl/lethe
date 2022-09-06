/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - 2019 by the Lethe authors
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

#ifndef lethe_gls_navier_stokes_h
#define lethe_gls_navier_stokes_h

#include <core/exceptions.h>

#include <solvers/copy_data.h>
#include <solvers/navier_stokes_base.h>
#include <solvers/navier_stokes_scratch_data.h>

#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>

using namespace dealii;

/**
 * @brief A solver class for the Navier-Stokes equation using GLS stabilization
 *
 * @tparam dim An integer that denotes the dimension of the space in which
 * the flow is solved
 *
 * @ingroup solvers
 */

template <int dim>
class GLSNavierStokesSolver
  : public NavierStokesBase<dim, TrilinosWrappers::MPI::Vector, IndexSet>
{
public:
  GLSNavierStokesSolver(SimulationParameters<dim> &nsparam);
  ~GLSNavierStokesSolver();

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
  setup_dofs_fd();

  /**
   * @brief update non zero constraint if the boundary is time dependant
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
   *  @brief Assembles the matrix associated with the solver without computing the preconditioner
   */
  virtual void
  assemble_system_matrix_without_preconditioner();

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
   * @brief Assemble the local matrix for a given cell.
   *
   * This function is used by the WorkStream class to assemble
   * the system matrix. It is a thread safe function.
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_matrix(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /**
   * @brief Assemble the local rhs for a given cell
   *
   * @param cell The cell for which the local matrix is assembled.
   *
   * @param scratch_data The scratch data which is used to store
   * the calculated finite element information at the gauss point.
   * See the documentation for NavierStokesScratchData for more
   * information
   *
   * @param copy_data The copy data which is used to store
   * the results of the assembly over a cell
   */
  virtual void
  assemble_local_system_rhs(
    const typename DoFHandler<dim>::active_cell_iterator &cell,
    NavierStokesScratchData<dim> &                        scratch_data,
    StabilizedMethodsTensorCopyData<dim> &                copy_data);

  /**
   * @brief sets up the vector of assembler functions
   */
  virtual void
  setup_assemblers();


  /**
   * @brief Copy local cell information to global matrix
   */

  virtual void
  copy_local_matrix_to_global_matrix(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);

  /**
   * @brief Copy local cell rhs information to global rhs
   */

  virtual void
  copy_local_rhs_to_global_rhs(
    const StabilizedMethodsTensorCopyData<dim> &copy_data);


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
                      const bool renewed_matrix = true);

private:
  /**
   * @brief Assembles an L2_projection matrix for the velocity and the pressure.
   * This L2 projection matrix can be used to set the initial condition for
   * the Navier-Stokes problem
   */
  void
  assemble_L2_projection();



  /**
   * @brief GMRES solver with ILU(N) preconditioning
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual);

  /**
   * @brief BiCGStab solver with ILU(N) preconditioning
   */
  void
  solve_system_BiCGStab(const bool   initial_step,
                        const double absolute_residual,
                        const double relative_residual);

  /**
   * @brief GMRES solver with AMG preconditioner with ILU smoother and coarsener
   */
  void
  solve_system_AMG(const bool   initial_step,
                   const double absolute_residual,
                   const double relative_residual);

  /**
   * @brief Direct solver using TrilinosWrappers::SolverDirect
   * The use of this solver should be avoided for 3D probelm
   */
  void
  solve_system_direct(const bool   initial_step,
                      const double absolute_residual,
                      const double relative_residual);

  /**
   * @brief  Set-up AMG preconditioner
   */
  void
  setup_AMG();

  /**
   * @brief Set-up ILU preconditioner
   */
  void
  setup_ILU();


  /**
   * Members
   */
protected:
  TrilinosWrappers::SparseMatrix system_matrix;

private:
  SparsityPattern                                    sparsity_pattern;
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_preconditioner;
  int current_preconditioner_fill_level;
  int initial_preconditioner_fill_level;
};


#endif
