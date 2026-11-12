// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom.h
 * @brief Auxiliary physics for the transport of a particle size distribution
 * with the quadrature method of moments (QMOM), coupled to the QMOMKit
 * library.
 *
 * The unknowns are two regular scalars, the ethanol concentration and the
 * free lipids concentration, and the first raw moments of the particle size
 * distribution. All of them are discretized together in a single FESystem of
 * discontinuous Galerkin elements of the same degree, so that every support
 * point carries the full local state.
 *
 * Each time step of the fluid dynamics is split into:
 *  - T: the physical-space transport of all the fields by the fluid velocity,
 *    solved by Lethe with an upwind discontinuous Galerkin discretization and
 *    an explicit strong-stability-preserving Runge-Kutta (SSP-RK) scheme,
 *    sub-cycled to satisfy the CFL condition;
 *  - P: the phase-space evolution (nucleation, growth and coalescence),
 *    applied independently at every support point through
 *    apply_qmom_sources(), which is delegated to QMOMKit.
 * The two are combined with Lie (T then P) or Strang (P/2, T, P/2)
 * splitting.
 *
 * The QMOM physics is advanced explicitly. It is therefore driven by the
 * linear solver strategy: assemble_system_matrix() and assemble_system_rhs()
 * do nothing, and solve_linear_system() advances the fields over the time
 * step.
 */

#ifndef lethe_qmom_h
#define lethe_qmom_h

#include <core/multiphysics.h>
#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/parameters_multiphysics.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/qmom_field_function.h>
#include <solvers/qmom_scratch_data.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/types.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
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
  InvalidNumberOfFluid,
  int,
  << "The QMOM physics is enabled, but the number of fluids is set to " << arg1
  << ". The QMOM physics solver only supports 1 fluid.");

/**
 * @brief Quadrature method of moments (QMOM), solved as an auxiliary physics
 * coupled to the fluid dynamics solver.
 *
 * Component layout of the FESystem: component ethanol_component holds the
 * ethanol concentration, component free_lipids_component holds the free
 * lipids concentration and components [first_moment_component,
 * n_components) hold the raw moments M_0 to M_{n_moments - 1} of the particle
 * size distribution.
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
   * @param[in] multiphysics_interface Multiphysics interface through which the
   * QMOM physics exchanges fields with the other physics.
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
   * @brief Gather the concentrations and the moments for output.
   *
   * @return Vector of output structures holding the QMOM solution.
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
   * @brief Shift the solution history.
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Post-process the QMOM solution: compute the L2 error of every
   * field if the analytical solution is enabled.
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
   * @brief Prepare the solution vectors for mesh refinement/coarsening.
   */
  void
  pre_mesh_adaptation() override;

  /**
   * @brief Interpolate the solution vectors onto the adapted mesh.
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
   * @brief Serialize the present and previous solutions for restart.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Deserialize the present and previous solutions for restart.
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
   * @brief Distribute the DoFs, allocate the vectors, and register the QMOM
   * fields in the multiphysics interface.
   */
  void
  setup_dofs() override;

  /**
   * @brief Interpolate the initial conditions of the "QMOM" subsection of the
   * "initial conditions".
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Advance the QMOM fields over the current time step. This is the
   * last function called by the linear solver strategy, the assembly
   * functions being empty for this explicit physics.
   */
  void
  solve_linear_system() override;

  /**
   * @brief Get the DoFHandler of the QMOM physics.
   *
   * @return Reference to the DoFHandler.
   */
  const DoFHandler<dim> &
  get_dof_handler() override
  {
    return *dof_handler;
  }

  /**
   * @brief Get the evaluation point.
   *
   * @return Reference to the evaluation point vector.
   */
  GlobalVectorType &
  get_evaluation_point() override
  {
    return evaluation_point;
  }

  /**
   * @brief Get the local (non-ghosted) evaluation point.
   *
   * @return Reference to the local evaluation point vector.
   */
  GlobalVectorType &
  get_local_evaluation_point() override
  {
    return local_evaluation_point;
  }

  /**
   * @brief Get the Newton update. Not used by this explicit physics.
   *
   * @return Reference to the Newton update vector.
   */
  GlobalVectorType &
  get_newton_update() override
  {
    return newton_update;
  }

  /**
   * @brief Get the present solution.
   *
   * @return Reference to the present solution vector.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    return *present_solution;
  }

  /**
   * @brief Get the system right-hand side. Not used by this explicit physics.
   *
   * @return Reference to the right-hand side vector.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    return system_rhs;
  }

  /**
   * @brief Get the non-zero constraints. They are empty, since the
   * discontinuous Galerkin discretization has no hanging-node constraints and
   * imposes the boundary conditions weakly.
   *
   * @return Reference to the non-zero constraints.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }

  /**
   * @brief Print the norms of the Newton update. Not used by this explicit
   * physics.
   */
  void
  output_newton_update_norms() override;

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
   * @brief Does nothing: the QMOM physics is advanced explicitly and has no
   * system matrix.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Does nothing: the transport operator is evaluated at every stage
   * of the time integration in solve_linear_system().
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Define the (empty) constraints of the discontinuous Galerkin
   * discretization.
   */
  void
  define_constraints();

  /**
   * @brief Advance the fields over the time step of the fluid dynamics with
   * the selected splitting of transport and phase-space evolution.
   */
  void
  advance_time_step();

  /**
   * @brief Advance the transport of the fields over a time step, with the
   * number of sub-steps required by the CFL condition.
   *
   * @param[in] time_step_start Time at the beginning of the step.
   * @param[in] time_step Duration of the step.
   */
  void
  advance_transport(const double time_step_start, const double time_step);

  /**
   * @brief Apply the phase-space evolution (nucleation, growth and
   * coalescence) at every locally owned support point.
   *
   * @param[in] time_step Duration of the phase-space step.
   */
  void
  apply_phase_space_step(const double time_step);

  /**
   * @brief Evaluate the transport operator L(u) = M^{-1} R(u) of the
   * semi-discrete system du/dt = L(u), where M is the block-diagonal mass
   * matrix and R the upwind discontinuous Galerkin residual of the advection
   * of all the fields.
   *
   * @tparam FluidVectorType Vector type of the fluid dynamics solution.
   * @param[in] solution Ghosted QMOM fields at which the operator is
   * evaluated.
   * @param[in] time Time at which the operator is evaluated (boundary
   * conditions and fluid velocity).
   * @param[in] interpolation_factor Position of time in the time step of the
   * fluid dynamics, between 0 and 1, used to interpolate the fluid velocity.
   * @param[in] previous_fluid_solution Fluid solution at the beginning of the
   * time step.
   * @param[in] present_fluid_solution Fluid solution at the end of the time
   * step.
   * @param[out] time_derivative Locally owned values of L(u).
   */
  template <typename FluidVectorType>
  void
  evaluate_transport_operator(const GlobalVectorType &solution,
                              const double            time,
                              const double            interpolation_factor,
                              const FluidVectorType  &previous_fluid_solution,
                              const FluidVectorType  &present_fluid_solution,
                              GlobalVectorType       &time_derivative);

  /**
   * @brief Assemble the upwind numerical flux of an interior face.
   *
   * @param[in] scratch_data Scratch data reinitialized on the face.
   * @param[out] face_rhs Contributions to the interface degrees of freedom.
   */
  void
  assemble_interior_face_flux(const QMOMScratchData<dim> &scratch_data,
                              Vector<double>             &face_rhs) const;

  /**
   * @brief Assemble the numerical flux of a boundary face according to its
   * boundary condition.
   *
   * @param[in,out] scratch_data Scratch data reinitialized on the face. Its
   * upwind_values buffer is overwritten.
   * @param[in] boundary_id Boundary id of the face.
   * @param[in,out] local_rhs Contributions to the degrees of freedom of the
   * cell.
   */
  void
  assemble_boundary_face_flux(QMOMScratchData<dim>    &scratch_data,
                              const types::boundary_id boundary_id,
                              Vector<double>          &local_rhs) const;

  /**
   * @brief Compute the L2 error of every field with respect to the analytical
   * solution at the current time.
   *
   * @return L2 errors, ordered as the components of the FESystem.
   */
  std::vector<double>
  calculate_L2_errors();

  /// Number of regular scalars
  const unsigned int n_scalars;

  /// Number of transported moments
  const unsigned int n_moments;

  /// Total number of components of the FESystem
  const unsigned int n_components;

  /// Component of the moment of order 0. The moment of order k is stored in
  /// component first_moment_component + k.
  const unsigned int first_moment_component;

  /// Multiphysics interface used to exchange fields with other physics
  MultiphysicsInterface<dim> *multiphysics;

  /// Timer of the QMOM physics
  TimerOutput computing_timer;

  /// Simulation parameters
  const SimulationParameters<dim> &simulation_parameters;

  /// Triangulation shared by all physics
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  /// Simulation control
  std::shared_ptr<SimulationControl> simulation_control;

  /// DoFHandler of the QMOM physics
  std::shared_ptr<DoFHandler<dim>> dof_handler;

  /// FESystem: n_scalars concentration components + n_moments moment
  /// components, all discretized with the same discontinuous element
  std::unique_ptr<FiniteElement<dim>> fe;

  /// Mapping, cell and face quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::unique_ptr<Quadrature<dim>>     cell_quadrature;
  std::unique_ptr<Quadrature<dim - 1>> face_quadrature;

  /// Properties of the aqueous and organic phases of the fluid, given to the
  /// phase-space step
  Parameters::QMOMLipidNanoparticleProduction lipid_nanoparticle_production;

  /// Names of the fields, ordered as the components of the FESystem
  std::vector<std::string> field_names;

  /// Initial conditions of all the fields
  std::shared_ptr<QMOMFieldFunction<dim>> initial_conditions;

  /// Values imposed by the Dirichlet boundary conditions, by boundary id
  std::map<types::boundary_id, std::shared_ptr<QMOMFieldFunction<dim>>>
    dirichlet_functions;

  /// Analytical solution of all the fields, used to compute the L2 errors
  std::shared_ptr<QMOMFieldFunction<dim>> analytical_solution;

  /// Locally owned and locally relevant DoFs
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  /// Solution and solver vectors. evaluation_point, local_evaluation_point,
  /// newton_update and system_rhs are required by the PhysicsSolver
  /// interface, but are not used by this explicit physics.
  GlobalVectorType                  evaluation_point;
  GlobalVectorType                  local_evaluation_point;
  GlobalVectorType                  newton_update;
  std::shared_ptr<GlobalVectorType> present_solution;
  GlobalVectorType                  system_rhs;
  AffineConstraints<double>         nonzero_constraints;

  /// Vectors of the explicit time integration. The ghosted stage solution is
  /// the state at which the transport operator is evaluated, the other
  /// vectors only hold locally owned values.
  GlobalVectorType ghosted_stage_solution;
  GlobalVectorType step_initial_solution;
  GlobalVectorType stage_solution;
  GlobalVectorType transport_residual;
  GlobalVectorType transport_time_derivative;

  /// Previous solutions, kept for checkpointing and for the other physics
  std::shared_ptr<std::vector<GlobalVectorType>> previous_solutions;

  /// Solution transfer for mesh adaptation and checkpointing
  std::unique_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;
  std::vector<SolutionTransfer<dim, GlobalVectorType>>
    previous_solutions_transfer;

  /// L2 errors of the fields with respect to the analytical solution
  ConvergenceTable error_table;
};

#endif
