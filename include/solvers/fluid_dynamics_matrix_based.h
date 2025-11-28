// SPDX-FileCopyrightText: Copyright (c) 2019-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_matrix_based_h
#define lethe_fluid_dynamics_matrix_based_h

#include <core/exceptions.h>
#include <core/vector.h>

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
class FluidDynamicsMatrixBased
  : public NavierStokesBase<dim, GlobalVectorType, IndexSet>
{
public:
  FluidDynamicsMatrixBased(SimulationParameters<dim> &nsparam);
  ~FluidDynamicsMatrixBased();

  /**
   * @brief solve Solves the Navier-Stokes problem
   *
   * This functions solves the problem defined by the Navier-Stokes parameter
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
   * @brief Update the sum_over_previous_stages variable to use it in the solver.
   *
   * @param stage An unsigned integer which gives the index of the current stage
   * @param method SDIRK method for now (BDF methods are single stage)
   */
  void
  multi_stage_preresolution(
    unsigned int                                      stage,
    Parameters::SimulationControl::TimeSteppingMethod method) override;

  /**
   * @brief Calculate \f$ k_i \f$ from \f$ u*_i : k_i = \frac{u*_i - u_n}{time_step*a_ii} - \sum{j=1}^{i-1} a_ij * k_j \f$
   * Update \f$ \sum_{i=1}^{N_{stages}} b_i * k_i \f$
   *
   * @param stage An unsigned integer which gives the index of the current stage
   * @param method SDIRK methods for now (BDF methods are single stage)
   */
  virtual void
  multi_stage_postresolution(
    unsigned int                                      stage,
    Parameters::SimulationControl::TimeSteppingMethod method,
    double                                            time_step) override;

  /**
   * @brief \f$ u_{n+1} = u_n + time_step * \sum_{i=1}^{N_{stages}} b_i * k_i \f$
   *
   * @param time_step The current time step to build the solution at the next time step
   */
  virtual void
  update_multi_stage_solution(double time_step) override;

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
   * overridden by the read_checkpoint functionality.
   *
   **/
  virtual void
  set_initial_condition_fd(
    Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    bool                                          restart = false) override;

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
  setup_preconditioner() override;

  /**
   * @brief Define the zero constraints used to solved the problem that change
   * with other physics' solutions.
   *
   * It differs from FluidDynamicsMatrixBased::define_zero_constraints as it
   * changes in time depending on the physics' solutions. Currently, it is only
   * used to constraint solid with the temperature's evolution.
   */
  void
  define_dynamic_zero_constraints();

  /**
   * @brief Update mortar configuration.
   *
   * When the rotor domain is rotated, the mortar cells need to be reinitialized
   * according to the new rotor-stator interface configuration.
   */
  virtual void
  update_mortar_configuration() override;

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
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data);

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
    NavierStokesScratchData<dim>                         &scratch_data,
    StabilizedMethodsTensorCopyData<dim>                 &copy_data);

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
   */
  void
  solve_linear_system() override;

private:
  /**
   * @brief Assembles an L2_projection matrix for the velocity and the pressure.
   * This L2 projection matrix can be used to set the initial condition for
   * the Navier-Stokes problem
   */
  void
  assemble_L2_projection();


  /**
   * @brief Solver for the L2 Projection linear system at the initial time step to apply this initial condition. It always uses the GMRES solver with ILU(N) by default.
   * 
   * @param[in] absolute_residual Absolute residual for the solver stopping criterion.
   * @param[in] relative_residual Relative residual for the solver stopping criterion.
   */

  void
  solve_L2_system(double absolute_residual, double relative_residual);


  /**
   * @brief GMRES solver with ILU(N) preconditioning or AMG preconditioning
   * @param[in] absolute_residual Absolute residual for the solver stopping criterion.
   * @param[in] relative_residual Relative residual for the solver stopping criterion.
   */
  void
  solve_system_GMRES(const double absolute_residual,
                     const double relative_residual);

  /**
   * @brief BiCGStab solver with ILU(N) preconditioning
   * @param[in] absolute_residual Absolute residual for the solver stopping criterion.
   * @param[in] relative_residual Relative residual for the solver stopping criterion.
   */
  void
  solve_system_BiCGStab(const double absolute_residual,
                        const double relative_residual);
  /**
   * @brief Direct solver using TrilinosWrappers::SolverDirect
   * The use of this solver should be avoided for 3D problems
   */
  void
  solve_system_direct(const double absolute_residual,
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
