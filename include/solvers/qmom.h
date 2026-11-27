// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom.h
 * @brief Auxiliary physics for the transport of a particle size distribution
 * with the quadrature method of moments (QMOM).
 *
 * The evolution in phase space is computed by the QMOMKit library.
 *
 * The unknowns are user-defined scalars and the moments of the particle
 * size distribution, called the components of the QMOM physics. All the
 * components are discretized with the same scalar discontinuous Galerkin
 * element, on a single DoFHandler, and each component is stored in its own
 * vector. Since all the components are advected by the same fluid velocity,
 * the matrix of their advection is the same: it is assembled and
 * preconditioned once per time step, then solved with one right-hand side per
 * component. The scalars may also diffuse, with the "scalar diffusivities" of
 * the QMOM subsection. The diffusion is discretized with the symmetric interior
 * penalty (SIPG) method, and the Dirichlet boundary conditions are imposed
 * weakly by the Nitsche method. The scalars that share a diffusivity share a
 * transport matrix (the advection matrix plus the diffusion matrix scaled by
 * the diffusivity) and its preconditioner.
 *
 * Each time step of the fluid dynamics, from t^n to t^{n+1}, is split (Lie
 * splitting) into the following steps, run by solve_linear_system() once the
 * fluid velocity at t^{n+1} is solved:
 *  1. the transport of the scalars in physical space by the fluid velocity,
 *     with their source terms, solved with an upwind discontinuous Galerkin
 *     discretization and an implicit BDF time integration
 *     (advance_scalars_physical_space());
 *  2. the transport of the moments in physical space, with the same
 *     discretization and their source terms (advance_moments_physical_space());
 *  3. the evolution in phase space, independently at every degree of freedom
 *     (advance_phase_space()), skipped when nucleation, growth and
 *     coalescence are all disabled:
 *     - the scalars keep their values at the end of step 1 during the whole
 *       phase-space step;
 *     - the moments at the end of step 2 are evolved by the nucleation, growth
 *       and coalescence, whose rates may depend on these scalars;
 *     - if the process couples the scalars to the particles, the scalars are
 *       then updated with the change of the moments (e.g. the free lipids
 *       consumed by the lipid nanoparticles).
 *
 * Steps 1 and 2 are independent: the right-hand sides of all the components
 * are assembled from their values at the beginning of the time step, so their
 * order does not change the solution. Both must precede step 3. The time step
 * ends with the phase-space evolution, so that the moments of every degree of
 * freedom with particles are rebuilt from a quadrature at the end of every
 * time step, and that the scalars include the mass exchanged with the
 * particles.
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
#include <solvers/copy_data.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/qmom_field_function.h>
#include <solvers/qmom_scratch_data.h>
#include <solvers/qmom_source_terms.h>
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
  unsigned int,
  << "The QMOM physics is enabled, but the number of fluids is set to " << arg1
  << ". The QMOM physics solver only supports 1 fluid.");

/**
 * @brief Exception raised when a PhysicsSolver function that requires a single
 * solution vector, or a Newton update, is called on the QMOM physics.
 */
DeclExceptionMsg(
  QMOMSingleVectorFunctionCalled,
  "The QMOM physics stores one solution vector per component and is solved "
  "with the linear solver strategy: this function should not be called.");

/**
 * @brief Quadrature method of moments (QMOM), solved as an auxiliary physics
 * coupled to the fluid dynamics solver.
 *
 * Component layout: component i < n_scalars is the scalar i and component
 * n_scalars + k is the moment M_k of the particle size distribution,
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
   * @brief Calculate the L2 error of every component, if an analytical
   * solution is given.
   *
   * The errors are written to a file and, if the analytical solution is
   * verbose, printed. The output of the solution itself is done by
   * attach_solution_to_output().
   *
   * @param[in] first_iteration Whether this is the first call, in which case
   * the errors are not calculated.
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
   *
   * @param[in] checkpoint_slot_label Label of the checkpoint slot.
   */
  void
  write_checkpoint(const std::string &checkpoint_slot_label) override;

  /**
   * @brief Deserialize the present and previous solutions of every
   * component.
   *
   * @param[in] checkpoint_slot_label Label of the checkpoint slot.
   */
  void
  read_checkpoint(const std::string &checkpoint_slot_label) override;

  /**
   * @brief Gather the post-processing tables to serialize. None yet.
   *
   * @return Empty vector of table handlers.
   */
  std::vector<OutputStructTableHandler>
  gather_tables() override;

  /**
   * @brief Distribute the DoFs and allocate the vectors and the transport
   * matrix.
   *
   * The DoFHandler and the mapping of the QMOM physics are registered in the
   * multiphysics interface, but not the solution vectors.
   */
  void
  setup_dofs() override;

  /**
   * @brief Interpolate the initial conditions of every component.
   *
   * They are given in the "QMOM" subsection of the "initial conditions" and
   * are interpolated at the support points of the element.
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Update the time of the Dirichlet boundary conditions if they are
   * time-dependent.
   *
   * The boundary conditions are imposed weakly, through the upwind fluxes, so
   * there are no constraints to rebuild.
   */
  void
  update_boundary_conditions() override;

  /**
   * @brief Build the preconditioner (ILU or AMG) of every transport matrix:
   * the advection matrix and, if the scalars diffuse, the matrix of each
   * distinct diffusivity.
   */
  void
  setup_preconditioner() override;

  /**
   * @brief Advance the components over the current time step with the Lie
   * splitting.
   *
   * The scalars and then the moments are transported in physical space. If
   * nucleation, growth or coalescence is enabled, the moments then evolve in
   * phase space with the scalars at the end of their transport, which ends
   * with the update of the scalars coupled to the particles (see the
   * description of the file). This is the last function
   * called by the linear solver strategy.
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
  }

  /**
   * @brief Get the constraints, which are empty.
   *
   * The discontinuous Galerkin discretization has no hanging-node constraints
   * and imposes the boundary conditions weakly.
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
    AssertThrow(false, QMOMSingleVectorFunctionCalled());
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
   * @brief Assemble the transport matrices of the components.
   *
   * The advection matrix contains the BDF time derivative, the upwind
   * discontinuous Galerkin advection and the outflow through the boundaries,
   * with the present velocity of the fluid dynamics. It transports the
   * moments and the scalars that do not diffuse. The matrix of each distinct
   * diffusivity of the scalars is the advection matrix plus the diffusion
   * matrix scaled by that diffusivity. The diffusion matrix is assembled at
   * the first call after setup_dofs().
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the right-hand side of every component.
   *
   * It contains the previous solutions of the BDF time derivative, the source
   * term and the inflow through the Dirichlet boundaries, and, for the
   * scalars that diffuse, the Nitsche terms of the Dirichlet boundaries.
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
   * @brief Assemble diffusion_matrix, the discontinuous Galerkin diffusion
   * operator with a unit diffusivity.
   *
   * It contains the cell term, the symmetric interior penalty (SIPG) terms of
   * the interior faces and the Nitsche terms of the Dirichlet boundaries. The
   * outlet and noflux boundaries have no diffusive flux. It does not depend on
   * the velocity or on the time step, and is only assembled when the degrees
   * of freedom change.
   */
  void
  assemble_diffusion_matrix();

  /**
   * @brief Add the contributions of a cell and of its interior faces to a
   * global matrix.
   *
   * A face between a locally owned and a ghost cell is assembled by the two
   * processes that own its cells. Each process only adds the rows that it
   * owns, so that every row receives the flux of each of its faces exactly
   * once.
   *
   * @param[in] copy_data Contributions of the cell and of its interior faces.
   * @param[in,out] matrix Global matrix to which the contributions are added.
   */
  void
  copy_local_to_global_matrix(const QMOMCopyData             &copy_data,
                              TrilinosWrappers::SparseMatrix &matrix) const;

  /**
   * @brief Define the (empty) constraints of the discontinuous Galerkin
   * discretization.
   */
  void
  define_constraints();

  /**
   * @brief Set up the ILU preconditioner of every transport matrix, stored in
   * ilu_preconditioners.
   */
  void
  setup_ilu();

  /**
   * @brief Set up the AMG preconditioner of every transport matrix, stored in
   * amg_preconditioners.
   */
  void
  setup_amg();

  /**
   * @brief Transport the moments in physical space over the time step (step 2
   * of the time step, see the description of the file).
   *
   * The transport system of every moment is solved, with the matrix and the
   * right-hand sides assembled by the linear solver strategy.
   */
  void
  advance_moments_physical_space();

  /**
   * @brief Transport the scalars in physical space over the time step (step 1
   * of the time step, see the description of the file).
   *
   * The transport system of every scalar is solved, with the right-hand sides
   * assembled by the linear solver strategy. Does nothing without scalars.
   *
   * The scalars that do not diffuse share the advection matrix and its
   * preconditioner with the moments. The others use the matrix of their
   * diffusivity, which includes the diffusion. Their source in physical space
   * is the source term of the "QMOM" subsection of the source term subsection:
   * the exchange of mass with the particles is applied at the end of the
   * phase-space step (see advance_phase_space()).
   */
  void
  advance_scalars_physical_space();

  /**
   * @brief Solve the transport system of one component over the time step.
   *
   * The system is solved with the transport matrix of the component and its
   * preconditioner, for the correction of the present solution of the
   * component: the right-hand side is the residual of that solution, and the
   * tolerance is relative to it.
   *
   * @param[in] component Component to transport: scalar component if
   * component < n_scalars, moment M_{component - n_scalars} otherwise.
   * @param[in] rescale_metric Metric used to rescale the residual, given by
   * get_residual_rescale_metric().
   */
  void
  solve_component_transport(const unsigned int component,
                            const double       rescale_metric);

  /**
   * @brief Evolve the moments in phase space at every locally owned degree of
   * freedom, with QMOMKit.
   *
   * This is step 3 of the time step (see the description of the file).
   *
   * It starts from the scalars and the moments at the end of their transport
   * in physical space.
   *
   * Since the element is nodal and discontinuous, the moments of a degree of
   * freedom are the moments at its support point: they are inverted into a
   * quadrature (Wheeler algorithm), the quadrature is evolved over the time
   * step by qmomkit::advance_qmom_phase_space(), and the moments are rebuilt
   * from it. The mechanisms are the callbacks of phase_space_source_terms:
   * the nodes move with the growth rate, the coalescence changes the moments
   * through the quadrature and the nucleation adds particles. They receive
   * the scalars of the degree of freedom, which are constant during the step.
   * Degrees of freedom without particles (zeroth moment below the cutoff of
   * QMOMKit) only change by nucleation. If the process couples the scalars
   * to the particles (QMOMSourceTerms::scalars_are_coupled()), the scalars
   * are then updated with the change of the moments (e.g. the free lipids
   * consumed by the lipid nanoparticles). If the moments of any degree of
   * freedom cannot be inverted or evolved, or if its particles consume more
   * of a scalar than available, an exception is thrown on every rank.
   *
   * @param[in] time_step Duration of the phase-space step, which ends at the
   * current time of the simulation control.
   */
  void
  advance_phase_space(const double time_step);

  /**
   * @brief Compute the L2 error of every component with respect to the
   * analytical solution at the current time, relative to the L2 norm of that
   * solution. The error of a component whose analytical solution vanishes is
   * absolute.
   *
   * @return Relative L2 errors, ordered as the components.
   */
  std::vector<double>
  calculate_L2_errors();

  /// Number of user-defined scalars
  const unsigned int n_scalars;

  /// Number of transported moments
  const unsigned int n_moments;

  /// Total number of components: n_scalars + n_moments
  const unsigned int n_components;

  /// Whether growth, nucleation, or coalescence is enabled. Otherwise, the
  /// phase-space step is skipped.
  const bool has_phase_space_mechanism;

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

  /// Growth, nucleation and coalescence of the particles in phase space
  std::unique_ptr<QMOMSourceTerms> phase_space_source_terms;

  /// Locally owned and locally relevant DoFs
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  /// Empty constraints of the discontinuous Galerkin discretization
  AffineConstraints<double> constraints;

  /// Diffusivity of each transport operator. The first operator, without
  /// diffusion, transports the moments and the scalars that do not diffuse.
  /// Each distinct positive diffusivity of the scalars has its own operator.
  std::vector<double> operator_diffusivities;

  /// Index of the transport operator of each component, in
  /// operator_diffusivities: component_operators[component]
  std::vector<unsigned int> component_operators;

  /// Transport matrix of each operator: system_matrices[0] is the advection
  /// matrix
  std::vector<TrilinosWrappers::SparseMatrix> system_matrices;

  /// Diffusion matrix with a unit diffusivity, only used if the scalars
  /// diffuse
  TrilinosWrappers::SparseMatrix diffusion_matrix;

  /// Whether diffusion_matrix is assembled for the present degrees of freedom
  bool diffusion_matrix_is_assembled = false;

  /// ILU preconditioner of each transport matrix
  std::vector<std::shared_ptr<TrilinosWrappers::PreconditionILU>>
    ilu_preconditioners;

  /// AMG preconditioner of each transport matrix
  std::vector<std::shared_ptr<TrilinosWrappers::PreconditionAMG>>
    amg_preconditioners;

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

  /// Relative L2 errors of the components with respect to the analytical
  /// solution
  ConvergenceTable error_table;
};

#endif
