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

template <typename VectorType>
class PreconditionBase;

/**
 * @brief A geometric multigrid preconditioner compatible with the
 * matrix-free solver.
 */
template <int dim>
class MFNavierStokesPreconditionGMG
{
  using VectorType     = LinearAlgebra::distributed::Vector<double>;
  using LSTransferType = MGTransferMatrixFree<dim, double>;
  using GCTransferType = MGTransferGlobalCoarsening<dim, VectorType>;
  using OperatorType   = NavierStokesOperatorBase<dim, double>;
  using SmootherPreconditionerType = PreconditionBase<VectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerTypeLS = PreconditionMG<dim, VectorType, LSTransferType>;
  using PreconditionerTypeGC = PreconditionMG<dim, VectorType, GCTransferType>;

public:
  /**
   * @brief Construct a new precondition GMG object. Sets constraints,
   * operators and transfer objects.
   *
   * @param[in] simulation_parameters Object containing all parameters specified
   * in input file.
   * @param[in] dof_handler Describes the layout of DoFs and the type of FE.
   * @param[in] dof_handler_fe_q_iso_q1 Describes the layout of DoFs for
   * FE_Q_iso_Q1 elements.
   * @param[in] mapping Describes the transformations from unit to real cell.
   * @param[in] cell_quadrature Required for local operations on cells.
   * @param[in] forcing_function Function specified in parameter file as source
   * term.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] fe Describes the FE system for the vector-valued problem.
   */
  MFNavierStokesPreconditionGMG(
    const SimulationParameters<dim>         &simulation_parameters,
    const DoFHandler<dim>                   &dof_handler,
    const DoFHandler<dim>                   &dof_handler_fe_q_iso_q1,
    const std::shared_ptr<Mapping<dim>>     &mapping,
    const std::shared_ptr<Quadrature<dim>>  &cell_quadrature,
    const std::shared_ptr<Function<dim>>     forcing_function,
    const std::shared_ptr<SimulationControl> simulation_control,
    const std::shared_ptr<FESystem<dim>>     fe);

  /**
   * @brief Initialize smoother, coarse grid solver and multigrid object
   * needed for the geometric multigrid preconditioner.
   *
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] flow_control Required for dynamic flow control.
   * @param[in] present_solution Previous solution needed to evaluate the non
   * linear term.
   * @param[in] time_derivative_previous_solutions Vector storing time
   * derivatives of previous solutions.
   */
  void
  initialize(const std::shared_ptr<SimulationControl> simulation_control,
             FlowControl<dim>                        &flow_control,
             const VectorType                        &present_solution,
             const VectorType &time_derivative_previous_solutions);

  /**
   * @brief Calls the v cycle function of the multigrid object.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Prints relevant multigrid information
   *
   */
  void
  print_relevant_info() const;

  /**
   * @brief Getter function for all level operators.
   *
   * @return Multigrid object that contains all level operators.
   */
  const MGLevelObject<std::shared_ptr<OperatorType>> &
  get_mg_operators() const;

private:
  /**
   * @brief Set the up AMG object needed for coarse-grid solver or
   * preconditioning.
   *
   */
  void
  setup_AMG();

  /**
   * @brief Set the up ILU object needed for coarse-grid solver or
   * preconditioning.
   *
   */
  void
  setup_ILU();

  /// Min level of the multigrid hierarchy
  unsigned int minlevel;

  /// Max level of the multigrid hierarchy
  unsigned int maxlevel;

  /// Intermediate level of the multigrid hierarchy
  unsigned int intlevel;

  /// Triangulations for the global coarsening case
  std::vector<std::shared_ptr<const Triangulation<dim>>>
    coarse_grid_triangulations;

  /// DoF handlers for each of the levels of the global coarsening algorithm
  MGLevelObject<DoFHandler<dim>> dof_handlers;

  /// Transfers for each of the levels of the global coarsening algorithm
  MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers;

  /// Level operators for the geometric multigrid
  MGLevelObject<std::shared_ptr<OperatorType>> mg_operators;

  /// Multigrid level object storing all operators
  std::shared_ptr<mg::Matrix<VectorType>> mg_matrix;

  /// Interface edge matrix needed only for local smoothing
  std::shared_ptr<mg::Matrix<VectorType>> mg_interface_matrix_in;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;

  /// Smoother object
  std::shared_ptr<
    MGSmootherPrecondition<OperatorType, SmootherType, VectorType>>
    mg_smoother;

  /// Collection of boundary constraints and refinement edge constrations for
  /// the different levels in the local smoothing approach.
  MGConstrainedDoFs mg_constrained_dofs;

  /// Transfer operator for local smoothing
  std::shared_ptr<LSTransferType> mg_transfer_ls;

  /// Transfer operator for global coarsening
  std::shared_ptr<GCTransferType> mg_transfer_gc;

  /// Algebraic multigrid as coarse grid solver
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> precondition_amg;

  /// Incomplete LU as coarse grid solver
  std::shared_ptr<TrilinosWrappers::PreconditionILU> precondition_ilu;

  /// Direct solver as coarse grid solver
  std::shared_ptr<TrilinosWrappers::SolverDirect> precondition_direct;

  /// Solver control for the coarse grid solver
  std::shared_ptr<ReductionControl> coarse_grid_solver_control;

  /// Solver control for the direct solver
  std::shared_ptr<SolverControl> direct_solver_control;

  /// GMRES as coarse grid solver
  std::shared_ptr<SolverGMRES<VectorType>> coarse_grid_solver;

  /// Multigrid wrapper for the coarse grid solver
  std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse;

  /// Solver control for the coarse grid solver (intermediate level)
  std::shared_ptr<SolverControl> coarse_grid_solver_control_intermediate;

  /// Multigrid wrapper for the coarse grid solver (intermediate level)
  std::shared_ptr<MGCoarseGridBase<VectorType>> mg_coarse_intermediate;

  /// GMRES as coarse grid solver (intermediate level)
  std::shared_ptr<SolverGMRES<VectorType>> coarse_grid_solver_intermediate;

  /// Multigrid method (intermediate level)
  std::shared_ptr<Multigrid<VectorType>> mg_intermediate;

  /// Global coarsening multigrid preconditioner object (intermediate level)
  std::shared_ptr<PreconditionMG<dim, VectorType, GCTransferType>>
    gc_multigrid_preconditioner_intermediate;

  /// Multigrid method
  std::shared_ptr<Multigrid<VectorType>> mg;

  /// Local smoothing multigrid preconditioner object
  std::shared_ptr<PreconditionMG<dim, VectorType, LSTransferType>>
    ls_multigrid_preconditioner;

  /// Global coarsening multigrid preconditioner object
  std::shared_ptr<PreconditionMG<dim, VectorType, GCTransferType>>
    gc_multigrid_preconditioner;

  /// Conditional Ostream
  ConditionalOStream pcout;

  /// Simulation parameters
  SimulationParameters<dim> simulation_parameters;

  /// Describes the layout of DoFs and the type of FE.
  const DoFHandler<dim> &dof_handler;

  /// Describes the layout of DoFs for FE_Q_iso_Q1 elements.
  const DoFHandler<dim> &dof_handler_fe_q_iso_q1;

  /// Vector holding number of coarse grid iterations
  mutable std::vector<unsigned int> coarse_grid_iterations;

public:
  /// Timer for specific geometric multigrid components.
  mutable TimerOutput mg_setup_timer;

  /// Internal timer for vmult timings
  mutable TimerOutput mg_vmult_timer;
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
  using VectorType = LinearAlgebra::distributed::Vector<double>;

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
   * @brief Set up appropriate preconditioner.
   *
   */
  void
  setup_preconditioner() override;

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
   */
  void
  setup_GMG();

  /**
   * @brief Setup the implicit LU preconditioner and call the solve function of the
   * linear solver. Attention: an actual matrix needs to be constructed using
   * the matrix-free operator.
   */
  void
  setup_ILU();

  /**
   * @brief Prints the setup times for the geometric multigrid preconditioner.
   *
   */
  void
  print_mg_setup_times();

  /**
   * @brief  Provide present and previous flow solutions to the multiphysics
   * interface. These solutions can be accessed through the multiphysics
   * interface by the different physics.
   */
  void
  update_solutions_for_multiphysics();

protected:
  /**
   * @brief Matrix-free operator in used for all the matrix-vector multiplications calls (vmult).
   *
   */
  std::shared_ptr<NavierStokesOperatorBase<dim, double>> system_operator;

  /**
   * @brief Geometric multigrid preconditioner.
   *
   */
  std::shared_ptr<MFNavierStokesPreconditionGMG<dim>> gmg_preconditioner;

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
   * @brief Describes the layout of DoFs for FE_Q_iso_Q1 elements.
   *
   */
  DoFHandler<dim> dof_handler_fe_q_iso_q1;

  /**
   * @brief Trilinos vector storing the average velocities that are provided to other
   * physics.
   *
   */
  TrilinosWrappers::MPI::Vector multiphysics_average_velocities;

  /**
   * @brief Trilinos vector storing the present solution that is provided to other physics.
   *
   */
  TrilinosWrappers::MPI::Vector multiphysics_present_solution;

  /**
   * @brief Vector storing trilinos vectors containing the previous solutions that are
   * provided to other physics.
   *
   */
  std::vector<TrilinosWrappers::MPI::Vector> multiphysics_previous_solutions;
};

#endif
