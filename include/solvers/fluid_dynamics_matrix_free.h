// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_matrix_free_h
#define lethe_fluid_dynamics_matrix_free_h

#include <core/exceptions.h>

#include <solvers/fluid_dynamics_matrix_free_operators.h>
#include <solvers/navier_stokes_base.h>

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

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

template <typename VectorType, typename VectorTypePrecondition>
class PreconditionAdapter;

namespace dealii
{
  template <class VectorType, class PreconditionerType>
  class MGCoarseGridApplyPreconditioner;
}

template <int dim, typename MGNumber>
class MGTransferMatrixFreeWrapper;

/**
 * @brief A geometric multigrid preconditioner compatible with the
 * matrix-free solver.
 */
template <int dim>
class MFNavierStokesPreconditionGMGBase
{
protected:
  using Number = double;

#ifndef LETHE_GMG_USE_FLOAT
  using MGNumber = double;
#else
#  if DEAL_II_VERSION_GTE(9, 7, 0)
  using MGNumber              = float;
#  else
  AssertThrow(
    false,
    ExcMessage(
      "Single precision for the geometric multigrid preconditioner requires a version of deal.II >= 9.7.0."));
#  endif
#endif

  using VectorType         = LinearAlgebra::distributed::Vector<Number>;
  using MGVectorType       = LinearAlgebra::distributed::Vector<MGNumber>;
  using TrilinosVectorType = LinearAlgebra::distributed::Vector<double>;
  using LSTransferType     = MGTransferMatrixFreeWrapper<dim, MGNumber>;
  using GCTransferType     = MGTransferGlobalCoarsening<dim, MGVectorType>;
  using OperatorType       = NavierStokesOperatorBase<dim, MGNumber>;
  using SmootherPreconditionerType = PreconditionBase<MGVectorType>;
  using SmootherType =
    PreconditionRelaxation<OperatorType, SmootherPreconditionerType>;
  using PreconditionerTypeLS =
    PreconditionMG<dim, MGVectorType, LSTransferType>;
  using PreconditionerTypeGC =
    PreconditionMG<dim, MGVectorType, GCTransferType>;

#if DEAL_II_VERSION_GTE(9, 7, 0)
  using CoarseGridSolverApply = MGCoarseGridApplyOperator<
    MGVectorType,
    PreconditionAdapter<MGVectorType, TrilinosVectorType>>;
#else
  using CoarseGridSolverApply = MGCoarseGridApplyPreconditioner<
    MGVectorType,
    PreconditionAdapter<MGVectorType, TrilinosVectorType>>;
#endif

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
   */
  MFNavierStokesPreconditionGMGBase(
    const SimulationParameters<dim> &simulation_parameters,
    const DoFHandler<dim>           &dof_handler,
    const DoFHandler<dim>           &dof_handler_fe_q_iso_q1);

  /**
   * @brief Default destructor.
   */

  virtual ~MFNavierStokesPreconditionGMGBase() = default;

  /**
   * @brief Initializes geometric multigrid preconditioner and pre-calculate terms that are constant.
   *
   * @param[in] mapping Describes the transformations from unit to real cell.
   * @param[in] cell_quadrature Required for local operations on cells.
   * @param[in] forcing_function Function specified in parameter file as source
   * term.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] physical_properties_manager Required to have the updated values
   * of physical properties in case of ramp or viscous initial conditions.
   * @param[in] fe Describes the FE system for the vector-valued problem.
   */
  void
  reinit(const std::shared_ptr<Mapping<dim>>      &mapping,
         const std::shared_ptr<Quadrature<dim>>   &cell_quadrature,
         const std::shared_ptr<Function<dim>>      forcing_function,
         const std::shared_ptr<SimulationControl> &simulation_control,
         const std::shared_ptr<PhysicalPropertiesManager>
                                             &physical_properties_manager,
         const std::shared_ptr<FESystem<dim>> fe);

  /**
   * @brief Creates the operator for a given level.
   *
   * @param[in] level Identifier of the level to be created
   */
  virtual void
  create_level_operator(const unsigned int level) = 0;

  /**
   * @brief Initialize smoother, coarse grid solver and multigrid object
   * needed for the geometric multigrid preconditioner.
   */
  void
  initialize();

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


  /**
   * @brief Getter function for all level smoother preconditioners.
   *
   * @return Multigrid object that contains all level smoother preconditioners.
   */
  const MGLevelObject<std::shared_ptr<PreconditionBase<MGVectorType>>> &
  get_mg_smoother_preconditioners() const;

protected:
  /**
   * @brief Set up AMG object needed for coarse-grid solver or
   * preconditioning.
   *
   */
  void
  setup_AMG();

  /**
   * @brief Set up ILU object needed for coarse-grid solver or
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
  MGLevelObject<MGTwoLevelTransfer<dim, MGVectorType>> transfers;

  /// Level operators for the geometric multigrid
  MGLevelObject<std::shared_ptr<OperatorType>> mg_operators;

  /// Multigrid level object storing all operators
  std::shared_ptr<mg::Matrix<MGVectorType>> mg_matrix;

  /// Interface edge matrix needed only for local smoothing
  std::shared_ptr<mg::Matrix<MGVectorType>> mg_interface_matrix_in;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_operators;
  MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<OperatorType>>
    ls_mg_interface_in;

  /// Preconditioners associated to smoothers
  mutable MGLevelObject<std::shared_ptr<PreconditionBase<MGVectorType>>>
    mg_smoother_preconditioners;

  /// Smoother object
  std::shared_ptr<
    MGSmootherPrecondition<OperatorType, SmootherType, MGVectorType>>
    mg_smoother;

  /// Collection of boundary constraints and refinement edge constrations for
  /// the different levels in the local smoothing approach.
  MGLevelObject<MGConstrainedDoFs> mg_constrained_dofs;

  /// Transfer operator for local smoothing
  std::shared_ptr<LSTransferType> mg_transfer_ls;

  /// Transfer operator for global coarsening
  std::shared_ptr<GCTransferType> mg_transfer_gc;

  /// Coarse grid solver (algebraic multigrid, ILU, direct solver)
  std::shared_ptr<PreconditionAdapter<MGVectorType, TrilinosVectorType>>
    coarse_grid_precondition;

  /// Solver control for the coarse grid solver
  std::shared_ptr<SolverControl> coarse_grid_solver_control;

  /// Solver control for the direct solver
  std::shared_ptr<SolverControl> direct_solver_control;

  /// GMRES as coarse grid solver
  std::shared_ptr<SolverGMRES<TrilinosVectorType>> coarse_grid_solver;

  /// Multigrid wrapper for the coarse grid solver
  std::shared_ptr<MGCoarseGridBase<MGVectorType>> mg_coarse;

  /// Solver control for the coarse grid solver (intermediate level)
  std::shared_ptr<SolverControl> coarse_grid_solver_control_intermediate;

  /// Multigrid wrapper for the coarse grid solver (intermediate level)
  std::shared_ptr<MGCoarseGridBase<MGVectorType>> mg_coarse_intermediate;

  /// GMRES as coarse grid solver (intermediate level)
  std::shared_ptr<SolverGMRES<MGVectorType>> coarse_grid_solver_intermediate;

  /// Multigrid method (intermediate level)
  std::shared_ptr<Multigrid<MGVectorType>> mg_intermediate;

  /// Global coarsening multigrid preconditioner object (intermediate level)
  std::shared_ptr<PreconditionMG<dim, MGVectorType, GCTransferType>>
    gc_multigrid_preconditioner_intermediate;

  /// Local smoothing multigrid preconditioner object (intermediate level)
  std::shared_ptr<PreconditionMG<dim, MGVectorType, LSTransferType>>
    ls_multigrid_preconditioner_intermediate;

  /// Multigrid method
  std::shared_ptr<Multigrid<MGVectorType>> mg;

  /// Local smoothing multigrid preconditioner object
  std::shared_ptr<PreconditionMG<dim, MGVectorType, LSTransferType>>
    ls_multigrid_preconditioner;

  /// Global coarsening multigrid preconditioner object
  std::shared_ptr<PreconditionMG<dim, MGVectorType, GCTransferType>>
    gc_multigrid_preconditioner;

  /// Conditional Ostream
  ConditionalOStream pcout;

  /// Simulation parameters
  SimulationParameters<dim> simulation_parameters;

  /// Describes the layout of DoFs and the type of FE.
  const DoFHandler<dim> &dof_handler;

  /// Describes the layout of DoFs for FE_Q_iso_Q1 elements.
  DoFHandler<dim> dof_handler_fe_q_iso_q1;

  /// Vector holding number of coarse grid iterations
  mutable std::vector<unsigned int> coarse_grid_iterations;

  /// Temperature DoF handlers for each of the levels of the global coarsening
  /// algorithm
  MGLevelObject<DoFHandler<dim>> temperature_dof_handlers;

  /// Temperature transfers for each of the levels of the global coarsening
  /// algorithm
  MGLevelObject<MGTwoLevelTransfer<dim, MGVectorType>> transfers_temperature;

  /// Transfer operator for global coarsening for the temperature
  std::shared_ptr<GCTransferType> mg_transfer_gc_temperature;

public:
  /// Timer for specific geometric multigrid components.
  mutable TimerOutput mg_setup_timer;

  /// Internal timer for vmult timings
  mutable TimerOutput mg_vmult_timer;
};

/**
 * @brief A geometric multigrid preconditioner implementation for
 * incompressible flow.
 */
template <int dim>
class MFNavierStokesPreconditionGMG
  : public MFNavierStokesPreconditionGMGBase<dim>
{
public:
  using VectorType =
    typename MFNavierStokesPreconditionGMGBase<dim>::VectorType;
  using MGVectorType =
    typename MFNavierStokesPreconditionGMGBase<dim>::MGVectorType;
  using MGNumber = typename MFNavierStokesPreconditionGMGBase<dim>::MGNumber;
  using GCTransferType = MGTransferGlobalCoarsening<dim, MGVectorType>;

  /**
   * Constructor.
   */
  MFNavierStokesPreconditionGMG(
    const SimulationParameters<dim> &simulation_parameters,
    const DoFHandler<dim>           &dof_handler,
    const DoFHandler<dim>           &dof_handler_fe_q_iso_q1);

  /*
   * @brief Creates the operator for one level
   *
   * @param[in] level Level to be created.
   */
  void
  create_level_operator(const unsigned int level) override;

  /*
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
  initialize(const std::shared_ptr<SimulationControl> &simulation_control,
             FlowControl<dim>                         &flow_control,
             const VectorType                         &present_solution,
             const VectorType &time_derivative_previous_solutions);

  /**
   * @brief Transfer the current temperature solution to the different multigrid levels and compute the buoyancy term in each of the multigrid operators.
   *
   * @param[in] temperature_dof_handler DoF Handler used for the heat transfer.
   * @param[in] temperature_present_solution Present solution of the temperature
   * as given by the multiphysics interface.
   */
  void
  initialize_auxiliary_physics(const DoFHandler<dim> &temperature_dof_handler,
                               const VectorType &temperature_present_solution);
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
class FluidDynamicsMatrixFree
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
  FluidDynamicsMatrixFree(SimulationParameters<dim> &nsparam);

  /**
   * @brief Destructor.
   *
   */
  ~FluidDynamicsMatrixFree();

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
   * @brief Set the initial conditions for the solver. If the simulation is
   * restarted from a checkpoint, the initial solution setting is bypassed
   * and the checkpoint is instead read.
   *
   * @param[in] initial_condition_type The type of initial condition to be
   * set.
   *
   * @param[in] restart A boolean that indicates if the simulation is being
   * restarted. If set to true, the initial conditions are never set, but are
   * instead overriden by the read_checkpoint functionality.
   *
   **/
  virtual void
  set_initial_condition_fd(
    Parameters::FluidDynamicsInitialConditionType initial_condition_type,
    bool                                          restart = false) override;

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
   * @brief Update mortar configuration.
   *
   * When the rotor domain is rotated, the mortar cells need to be reinitialized
   * according to the new rotor-stator interface configuration.
   */
  virtual void
  update_mortar_configuration() override;

  void
  reinit_mortar_operators_mf();

  /**
   * @brief  Update the average velocity field solution in the multiphyscics interface.
   */
  virtual void
  update_multiphysics_time_average_solution() override;


  /**
   * @brief  Provide present and previous flow solutions to the multiphysics
   * interface. These solutions can be accessed through the multiphysics
   * interface by the different physics.
   */
  virtual void
  update_solutions_for_multiphysics() override;


  /**
   * @brief  Provide present multiphysics solutions from the multiphysics
   * interface to the fluid dynamics.
   */
  virtual void
  update_solutions_for_fluid_dynamics() override;

  /**
   * @brief Calculate and store time derivatives of previous solutions according to
   * time-stepping scheme to use them in the operator.
   */
  void
  calculate_time_derivative_previous_solutions();

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
   * linear system. If this is the case, the non_zero version of the
   * constraints are used for the Dirichlet boundary conditions.
   *
   * @param[in] renewed_matrix Indicates if the matrix has been reassembled,
   * and thus the preconditioner needs to be reassembled.
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
   * linear system. If this is the case, the non_zero version of the
   * constraints are used for the Dirichlet boundary conditions
   *
   * @param[in] absolute_residual Used to define the linear solver tolerance.
   *
   * @param[in] relative_residual Used to define the linear solver tolerance.
   */
  void
  solve_system_GMRES(const bool   initial_step,
                     const double absolute_residual,
                     const double relative_residual);

  void
  solve_system_direct(const bool   initial_step,
                      const double absolute_residual,
                      const double relative_residual);
  /**
   * @brief  Create the geometric multigrid preconditioner.
   */
  virtual void
  create_GMG();

  /**
   * @brief  Intialize the geometric multigrid preconditioner.
   */
  virtual void
  initialize_GMG();

  /**
   * @brief  Setup the geometric multigrid preconditioner.
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



protected:
  /**
   * @brief Prints the setup times for the geometric multigrid preconditioner.
   *
   */
  void
  print_mg_setup_times();

  /**
   * @brief Matrix-free operator in used for all the matrix-vector multiplications calls (vmult).
   *
   */
  std::shared_ptr<NavierStokesOperatorBase<dim, double>> system_operator;

  /**
   * @brief Geometric multigrid preconditioner.
   *
   */
  std::shared_ptr<MFNavierStokesPreconditionGMGBase<dim>> gmg_preconditioner;

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

  /**
   * @brief Vector storing the present temperature solution from the heat transfer solver.
   *
   */
  VectorType temperature_present_solution;

  /**
   * @brief Pointer to the physical properties object
   *
   */
  std::shared_ptr<PhysicalPropertiesManager> physical_properties_manager;
};

#endif
