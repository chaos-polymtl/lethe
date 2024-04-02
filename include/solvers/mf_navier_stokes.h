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

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>


using namespace dealii;


template <int dim>
class MFNavierStokesPreconditionGMG
{
  using VectorType     = LinearAlgebra::distributed::Vector<double>;
  using LSTransferType = MGTransferMatrixFree<dim, double>;
  using GCTransferType = MGTransferGlobalCoarsening<dim, VectorType>;
  using OperatorType   = NavierStokesOperatorBase<dim, double>;
  using SmootherPreconditionerType = DiagonalMatrix<VectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerTypeLS = PreconditionMG<dim, VectorType, LSTransferType>;
  using PreconditionerTypeGC = PreconditionMG<dim, VectorType, GCTransferType>;

public:
  void
  initialize_ls(
    TimerOutput                             &computing_timer,
    const DoFHandler<dim>                   &dof_handler,
    const SimulationParameters<dim>         &simulation_parameters,
    const std::shared_ptr<Mapping<dim>>     &mapping,
    const std::shared_ptr<FESystem<dim>>     fe,
    TimerOutput                             &mg_computing_timer,
    const std::shared_ptr<Quadrature<dim>>  &cell_quadrature,
    const std::shared_ptr<Function<dim>>     forcing_function,
    const VectorType                        &present_solution,
    const VectorType                        &time_derivative_previous_solutions,
    const ConditionalOStream                &pcout,
    const std::shared_ptr<SimulationControl> simulation_control) const;

  void
  initialize_gc(
    TimerOutput                             &computing_timer,
    const DoFHandler<dim>                   &dof_handler,
    const SimulationParameters<dim>         &simulation_parameters,
    const std::shared_ptr<Mapping<dim>>     &mapping,
    const std::shared_ptr<FESystem<dim>>     fe,
    TimerOutput                             &mg_computing_timer,
    const std::shared_ptr<Quadrature<dim>>  &cell_quadrature,
    const std::shared_ptr<Function<dim>>     forcing_function,
    const VectorType                        &present_solution,
    const VectorType                        &time_derivative_previous_solutions,
    const ConditionalOStream                &pcout,
    const std::shared_ptr<SimulationControl> simulation_control) const;

  void
  vmult(VectorType &dst, const VectorType &src) const;

  void
  clear() const;

private:
  // onlt GC
  mutable MGLevelObject<DoFHandler<dim>>                     dof_handlers;
  mutable MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;

  // level matrices
  mutable MGLevelObject<std::shared_ptr<OperatorType>>
                                                  mg_operators; // TODO: reuse
  mutable std::shared_ptr<mg::Matrix<VectorType>> mg_matrix;

  // edge matrices (only LS)
  mutable std::shared_ptr<mg::Matrix<VectorType>> mg_interface_matrix_in;
  mutable std::shared_ptr<mg::Matrix<VectorType>> mg_interface_matrix_out;
  mutable MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  mutable MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;
  mutable MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
                                                       ls_mg_interface_out;
  mutable MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_in;
  mutable MGLevelObject<std::shared_ptr<OperatorType>> mg_interface_out;

  // smoother
  mutable std::shared_ptr<
    MGSmootherPrecondition<OperatorType, SmootherType, VectorType>>
    mg_smoother;

  // transfer operators
  mutable MGConstrainedDoFs               mg_constrained_dofs;
  mutable std::shared_ptr<LSTransferType> mg_transfer_ls; // TODO: reuse
  mutable std::shared_ptr<GCTransferType> mg_transfer_gc; // TODO: reuse

  // coarse-grid solvers
  mutable TrilinosWrappers::PreconditionAMG        precondition_amg;
  mutable TrilinosWrappers::PreconditionILU        precondition_ilu;
  mutable std::shared_ptr<ReductionControl>        coarse_grid_solver_control;
  mutable std::shared_ptr<SolverGMRES<VectorType>> coarse_grid_solver;
  mutable std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  // multigrid as precondtioner
  mutable std::shared_ptr<Multigrid<VectorType>> mg;
  mutable std::shared_ptr<PreconditionMG<dim, VectorType, LSTransferType>>
    ls_multigrid_preconditioner;
  mutable std::shared_ptr<PreconditionMG<dim, VectorType, GCTransferType>>
    gc_multigrid_preconditioner;
};


/**
 * @brief A solver for the incompressible Navier-Stokes equations implemented
 * in a matrix-free fashion.
 *
 * A matrix-free stabilized solver for the incompressible Navier-Stokes
 * equations implemented in a matrix-free fashion with a sum-factorization
 * approach. It uses a continuous Galerkin discretization and solves the
 * fully-coupled discretized problem in a monolithic way using multigrid
 * preconditioners.
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
  using VectorType     = LinearAlgebra::distributed::Vector<double>;
  using LSTransferType = MGTransferMatrixFree<dim, double>;
  using GCTransferType = MGTransferGlobalCoarsening<dim, VectorType>;

public:
  /**
   * @brief Constructor that sets the finite element degree and system operator
   * according to simulation parameters.
   *
   * @param[in] nsparam Relevant parameters for the solver.
   */
  MFNavierStokesSolver(SimulationParameters<dim> &nsparam);

  /**
   * @brief Destructor.
   *
   */
  ~MFNavierStokesSolver();

  /**
   * @brief Solve the problem defined by simulation parameters by iterating
   * through time or through the mesh refinements.
   */
  virtual void
  solve();

protected:
  /**
   * @brief Setup the degrees of freedom, system constraints, system operator
   * and solution vectors.
   */
  virtual void
  setup_dofs_fd() override;

  /**
   * @brief Update non zero constraints if the boundary is time dependent.
   */
  void
  update_boundary_conditions();

  /**
   * @brief Set the initial conditions for the solver. If the simulation is
   * restarted from a checkpoint, the initial solution setting is bypassed
   * and the checkpoint is instead read.
   *
   * @param[in] initial_condition_type The type of initial condition to be set.
   *
   * @param[in] restart A boolean that indicates if the simulation is being
   * restarted. If set to true, the initial conditions are never set, but are
   * instead overriden by the read_checkpoint functionality.
   *
   **/
  virtual void
  set_initial_condition_fd(
    Parameters::InitialConditionType initial_condition_type,
    bool                             restart = false) override;

  /**
   * @brief Assemble the matrix associated with the solver. Only required for
   * compilation and it is not used for the matrix free solver.
   */
  virtual void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the system right hand side associated with the solver.
   */
  virtual void
  assemble_system_rhs() override;

  /**
   * @brief  Update the average velocity field solution in the multiphyscics interface.
   * Currently not implemented for this solver but required for compilation.
   */
  virtual void
  update_multiphysics_time_average_solution() override;

  /**
   * @brief Calculate and store time derivatives of previous solutions according to
   * time-stepping scheme to use them in the operator.
   */
  void
  calculate_time_derivative_previous_solutions();

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
   * @brief Solve the linear system of equations using the method specified in
   * the simulation parameters.
   *
   * @param[in] initial_step Indicates if this is the first solution of the
   * linear system. If this is the case, the non_zero version of the constraints
   * are used for the Dirichlet boundary conditions.
   *
   * @param[in] renewed_matrix Indicates if the matrix has been reassembled, and
   * thus the preconditioner needs to be reassembled.
   */
  void
  solve_linear_system(const bool initial_step,
                      const bool renewed_matrix = true) override;

private:
  /**
   * @brief Assemble a L2 projection matrix for the velocity and the pressure,
   * which can be used to set the initial condition for the Navier-Stokes
   * problem. Currently not implemented for this solver.
   */
  void
  assemble_L2_projection();

  /**
   * @brief GMRES solver with preconditioning.
   *
   * @param[in] initial_step Indicates if this is the first solution of the
   * linear system. If this is the case, the non_zero version of the constraints
   * are used for the Dirichlet boundary conditions
   *
   * @param[in] absolute_residual Used to define the linear solver tolerance.
   *
   * @param[in] relative_residual Used to define the linear solver tolerance.
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual);

  /**
   * @brief  Setup the geometric multigrid preconditioner and call the solve
   * function of the linear solver.
   *
   * @param[in] solver Linear solver object that needs the multigrid
   * preconditioner.
   */
  void
  solve_with_GMG(SolverGMRES<VectorType> &solver);

  /**
   * @brief Setup the implicit LU preconditioner and call the solve function of the
   * linear solver. Attention: an actual matrix needs to be constructed using
   * the matrix-free operator.
   */
  void
  solve_with_ILU(SolverGMRES<VectorType> &solver);

protected:
  /**
   * @brief Matrix-free operator in used for all the matrix-vector multiplications calls (vmult).
   *
   */
  std::shared_ptr<NavierStokesOperatorBase<dim, double>> system_operator;

  /**
   * @brief Geometric local smoothing multigrid preconditioner.
   *
   */
  std::shared_ptr<PreconditionMG<dim, VectorType, LSTransferType>>
    ls_multigrid_preconditioner;

  /**
   * @brief Geometric global coarsening multigrid preconditioner.
   *
   */
  std::shared_ptr<PreconditionMG<dim, VectorType, GCTransferType>>
    gc_multigrid_preconditioner;

  /**
   * @brief Implicit LU preconditioner.
   *
   */
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /**
   * @brief Vector to store the time derivatives of the previous solutions at the end
   * of each time step for time-dependent simulations.
   *
   */
  VectorType time_derivative_previous_solutions;

  /**
   * @brief Timer for specific geometric multigrid components.
   *
   */
  TimerOutput mg_computing_timer;
};

#endif
