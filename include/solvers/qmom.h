// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom.h
 * @brief Auxiliary physics for the transport of a particle size distribution
 * with the quadrature method of moments (QMOM), coupled to the QMOMKit
 * library.
 *
 * The unknowns are regular scalars and the first raw moments of the particle
 * size distribution, called the components of the QMOM physics. All the
 * components are discretized with the same scalar discontinuous Galerkin
 * element, on a single DoFHandler, and each component is stored in its own
 * vector. Since all the components are advected by the same fluid velocity,
 * the matrix of their transport is the same: it is assembled and
 * preconditioned once per time step, then solved with one right-hand side per
 * component.
 *
 * Each time step of the fluid dynamics is split (Lie splitting) into:
 *  - the transport of all the components in physical space by the fluid
 *    velocity, with their source terms, solved with an upwind discontinuous
 *    Galerkin discretization and an implicit BDF time integration
 *    (advance_physical_space());
 *  - the evolution of the moments in phase space (nucleation, growth and
 *    coalescence), applied independently at every degree of freedom
 *    (advance_phase_space()).
 *
 * The transport is linear in the components, since the fluid velocity is
 * given. The QMOM physics is therefore driven by the linear solver strategy:
 * assemble_system_matrix() assembles the transport matrix,
 * setup_preconditioner() builds its preconditioner, assemble_system_rhs()
 * assembles the right-hand side of every component and solve_linear_system()
 * advances the components over the time step.
 */

#ifndef lethe_qmom_h
#define lethe_qmom_h

#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/qmom_field_function.h>
#include <solvers/qmom_scratch_data.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/solution_transfer.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>


DeclException1(QMOMDimensionalityError,
               std::string,
               << "The QMOM physics only supports SI units, but the "
               << "dimensionality subsection sets " << arg1
               << ". Remove the dimensionality subsection or set these "
               << "parameters to 1.");

DeclException1(
  QMOMInvalidNumberOfFluids,
  int,
  << "The QMOM physics is enabled, but the number of fluids is set to " << arg1
  << ". The QMOM physics solver only supports 1 fluid.");

/**
 * @brief Quadrature method of moments (QMOM), solved as an auxiliary physics
 * coupled to the fluid dynamics solver.
 *
 * Component layout: component i < n_scalars is the scalar i and component
 * n_scalars + k is the raw moment M_k of the particle size distribution,
 * k < n_moments. This is the ordering of the components of the
 * QMOMFieldFunction of the initial conditions, boundary conditions, source
 * term and analytical solution.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class QMOM : public AuxiliaryPhysics<dim, GlobalVectorType>
{
public:
  /**
   * @brief Constructor of the QMOM physics. It verifies that the QMOM
   * parameters are consistent with the features supported so far.
   *
   * @param[in] multiphysics_interface Multiphysics interface, through which
   * the QMOM physics gets the fluid velocity.
   * @param[in] p_simulation_parameters Simulation parameters.
   * @param[in] p_triangulation Distributed triangulation shared by all physics.
   * @param[in] p_simulation_control Simulation control (time stepping).
   */
  QMOM(MultiphysicsInterface<dim>      *multiphysics_interface,
       const SimulationParameters<dim> &p_simulation_parameters,
       const std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                                &p_triangulation,
       const std::shared_ptr<SimulationControl> &p_simulation_control);

  /**
   * @brief Gather the components for output, one output structure per
   * component, since each component is stored in its own vector.
   *
   * @return Vector of output structures holding the QMOM components.
   */
  std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;

  /**
   * @brief Finish the simulation: print the L2 error table, if the
   * analytical solution is enabled, and the timers.
   */
  void
  finish_simulation() override;

  /**
   * @brief Shift the solution history of every component.
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Post-process the QMOM components: compute the L2 error of every
   * component if the analytical solution is enabled.
   *
   * @param[in] first_iteration Whether this is the first iteration.
   */
  void
  postprocess(bool first_iteration) override;

  /**
   * @brief Print the timer summary at each iteration.
   */
  void
  output_per_iteration_timer() override;

  /**
   * @brief Prepare the present and previous solutions of every component for
   * mesh refinement/coarsening.
   */
  void
  pre_mesh_adaptation() override;

  /**
   * @brief Interpolate the present and previous solutions of every component
   * onto the adapted mesh.
   */
  void
  post_mesh_adaptation() override;

  /**
   * @brief Compute the error estimate used for mesh adaptation. Not supported
   * yet for the QMOM physics.
   *
   * @param[in] ivar Variable and adaptation parameters.
   * @param[out] estimated_error_per_cell Estimated error per cell.
   */
  void
  compute_error_estimate(
    const std::pair<const Variable, Parameters::MultipleAdaptationParameters>
                          &ivar,
    dealii::Vector<float> &estimated_error_per_cell) override;

  /**
   * @brief Serialize the present and previous solutions of every component
   * for restart.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Deserialize the present and previous solutions of every
   * component.
   */
  void
  read_checkpoint() override;

  /**
   * @brief Gather the post-processing tables to serialize. None yet.
   *
   * @return Empty vector of table handlers.
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Distribute the DoFs, allocate the vectors and the transport
   * matrix, and register the DoFHandler and the mapping of the QMOM physics
   * in the multiphysics interface. The solution vectors are not registered.
   */
  void
  setup_dofs() override;

  /**
   * @brief Interpolate the initial conditions of the "QMOM" subsection of the
   * "initial conditions" at the support points of the element, component by
   * component.
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Update the time of the Dirichlet boundary conditions if they are
   * time-dependent. The boundary conditions are imposed weakly, through the
   * upwind fluxes, so there are no constraints to rebuild.
   */
  void
  update_boundary_conditions() override;

  /**
   * @brief Build the preconditioner (ILU or AMG) of the transport matrix,
   * shared by all the components.
   */
  void
  setup_preconditioner() override;

  /**
   * @brief Advance the components over the current time step with Lie
   * splitting:
   * transport in physical space, then evolution in phase space. This is the
   * last function called by the linear solver strategy.
   */
  void
  solve_linear_system() override;

  /**
   * @brief Get the DoFHandler shared by all the QMOM components.
   *
   * @return Reference to the DoFHandler.
   */
  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return *dof_handler;
  }

  /**
   * @brief Not available: the QMOM physics is solved with the linear solver
   * strategy and has no evaluation point.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Not available: the QMOM physics is solved with the linear solver
   * strategy and has no evaluation point.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_local_evaluation_point() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Not available: the QMOM physics is solved with the linear solver
   * strategy and has no Newton update.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_newton_update() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Not available: the QMOM physics stores one solution vector per
   * component.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Not available: the QMOM physics has one right-hand side per
   * component.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Get the constraints. They are empty, since the discontinuous
   * Galerkin discretization has no hanging-node constraints and imposes the
   * boundary conditions weakly.
   *
   * @return Reference to the constraints.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return constraints;
  }

  /**
   * @brief Not available: the QMOM physics is solved with the linear solver
   * strategy and has no Newton update.
   */
  void
  output_newton_update_norms() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Get the metric used to rescale the residual.
   *
   * @return Square root of the domain volume if residual rescaling is
   * enabled, 1 otherwise.
   */
  double
  get_residual_rescale_metric() const override;

private:
  /**
   * @brief Assemble the transport matrix shared by all the components: BDF
   * time derivative, upwind discontinuous Galerkin advection and outflow
   * through the boundaries, with the present velocity of the fluid dynamics.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the right-hand side of every component: previous
   * solutions of the BDF time derivative, source term and inflow through the
   * Dirichlet boundaries.
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Gather the present fluid velocity at the quadrature points of a
   * cell, from the block or non-block fluid dynamics solution.
   *
   * @param[in] velocity_cell Cell of the fluid dynamics DoFHandler.
   * @param[in,out] scratch_data Scratch data in which the velocity is stored.
   */
  void
  reinit_velocity_with_adequate_solution(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    QMOMScratchData<dim>                                 &scratch_data);

  /**
   * @brief Gather the present fluid velocity at the quadrature points of a
   * face, from the block or non-block fluid dynamics solution.
   *
   * @param[in] velocity_cell Cell of the fluid dynamics DoFHandler.
   * @param[in] face_no Face number in velocity_cell.
   * @param[in] sub_face_no Subface number in velocity_cell, or
   * numbers::invalid_unsigned_int if the face is not refined from this side.
   * @param[in,out] scratch_data Scratch data in which the velocity is stored.
   */
  void
  reinit_face_velocity_with_adequate_solution(
    const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
    const unsigned int                                    face_no,
    const unsigned int                                    sub_face_no,
    QMOMScratchData<dim>                                 &scratch_data);

  /**
   * @brief Define the (empty) constraints of the discontinuous Galerkin
   * discretization.
   */
  void
  define_constraints();

  /**
   * @brief Set up the ILU preconditioner stored in ilu_preconditioner.
   */
  void
  setup_ilu();

  /**
   * @brief Set up the AMG preconditioner stored in amg_preconditioner.
   */
  void
  setup_amg();

  /**
   * @brief Transport all the components in physical space over the time
   * step: solve the transport system of every component, with the matrix and
   * the right-hand sides assembled by the linear solver strategy.
   */
  void
  advance_physical_space();

  /**
   * @brief Evolve the moments in phase space (nucleation, growth and
   * coalescence) at every locally owned degree of freedom. Not implemented
   * yet: the moments are constant in phase space.
   *
   * @param[in] time_step Duration of the phase-space step.
   */
  void
  advance_phase_space(const double time_step);

  /**
   * @brief Compute the L2 error of every component with respect to the
   * analytical solution at the current time.
   *
   * @return L2 errors, ordered as the components.
   */
  std::vector<double>
  calculate_L2_errors();

  /// Error message of the PhysicsSolver functions that require a single
  /// solution vector
  static constexpr const char *single_vector_error_message =
    "The QMOM physics stores one solution vector per component and is solved "
    "with the linear solver strategy: this function should not be called.";

  /// Number of regular scalars
  const unsigned int n_scalars;

  /// Number of transported moments
  const unsigned int n_moments;

  /// Total number of components: n_scalars + n_moments
  const unsigned int n_components;

  /// Multiphysics interface, used to get the fluid velocity
  MultiphysicsInterface<dim> *multiphysics;

  /// Timer of the QMOM physics
  TimerOutput computing_timer;

  /// Simulation parameters
  const SimulationParameters<dim> &simulation_parameters;

  /// Triangulation shared by all physics
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  /// Simulation control
  std::shared_ptr<SimulationControl> simulation_control;

  /// DoFHandler shared by all the components
  std::shared_ptr<DoFHandler<dim>> dof_handler;

  /// Scalar discontinuous element shared by all the components
  std::unique_ptr<FiniteElement<dim>> fe;

  /// Mapping, cell and face quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::unique_ptr<Quadrature<dim>>     cell_quadrature;
  std::unique_ptr<Quadrature<dim - 1>> face_quadrature;

  /// Names of the components, used for the output and the error table
  std::vector<std::string> component_names;

  /// Initial conditions of all the components
  std::shared_ptr<QMOMFieldFunction<dim>> initial_conditions;

  /// Values imposed by the Dirichlet boundary conditions, by boundary id
  std::map<types::boundary_id, std::shared_ptr<QMOMFieldFunction<dim>>>
    dirichlet_functions;

  /// Source term of all the components
  std::shared_ptr<QMOMFieldFunction<dim>> source_term;

  /// Analytical solution of all the components, used to compute the L2
  /// errors
  std::shared_ptr<QMOMFieldFunction<dim>> analytical_solution;

  /// Locally owned and locally relevant DoFs
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  /// Empty constraints of the discontinuous Galerkin discretization
  AffineConstraints<double> constraints;

  /// Transport matrix, shared by all the components
  TrilinosWrappers::SparseMatrix system_matrix;

  /// ILU preconditioner of the transport matrix
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /// AMG preconditioner of the transport matrix
  std::shared_ptr<TrilinosWrappers::PreconditionAMG> amg_preconditioner;

  /// Ghosted present solution of each component:
  /// present_solutions[component]
  std::vector<GlobalVectorType> present_solutions;

  /// Ghosted previous solutions of each component, used by the BDF schemes:
  /// previous_solutions[component][p] is the solution p + 1 time steps ago
  std::vector<std::vector<GlobalVectorType>> previous_solutions;

  /// Right-hand side of each component: system_rhs[component]
  std::vector<GlobalVectorType> system_rhs;

  /// Solution transfer of the present and previous solutions of all the
  /// components, for mesh adaptation and checkpointing
  std::shared_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;

  /// L2 errors of the components with respect to the analytical solution
  ConvergenceTable error_table;
};

#endif
