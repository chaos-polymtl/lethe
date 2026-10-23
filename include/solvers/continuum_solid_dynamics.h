// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file continuum_solid_dynamics.h
 * @brief Auxiliary physics for the solid phase of the Euler-Euler (two-fluid)
 * model.
 *
 * The solid phase is treated as an interpenetrating continuum. Its unknowns
 * are the solid velocity (dim components) and the solid volume fraction
 * (1 component), discretized together in a single FESystem. The fluid phase
 * is solved by the existing volume-averaged Navier-Stokes (VANS) solver.
 *
 * In this first version, the class only provides the infrastructure (degrees
 * of freedom, solution vectors, output, mesh adaptation, checkpointing and
 * registration in the MultiphysicsInterface). The governing equations are not
 * assembled yet, so should_solve_auxiliary_physics() returns false.
 */

#ifndef lethe_continuum_solid_dynamics_h
#define lethe_continuum_solid_dynamics_h

#include <core/multiphysics.h>
#include <core/output_struct.h>
#include <core/parameters.h>
#include <core/simulation_control.h>
#include <core/vector.h>

#include <solvers/auxiliary_physics.h>
#include <solvers/multiphysics_interface.h>
#include <solvers/simulation_parameters.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/solution_transfer.h>

#include <memory>
#include <vector>

/**
 * @brief Exception thrown when the solid phase is enabled with a fluid
 * dynamics solver that does not solve the volume-averaged Navier-Stokes
 * equations.
 */
DeclExceptionMsg(
  ContinuumSolidDynamicsWithoutVANSError,
  "The \"continuum solid dynamics\" physics (Euler-Euler model) is enabled, but the fluid "
  "dynamics solver does not solve the volume-averaged Navier-Stokes (VANS) "
  "equations. The solid phase can only be used with a VANS solver.");

/**
 * @brief Solid phase of the Euler-Euler model, solved as an auxiliary physics
 * coupled to the VANS fluid dynamics solver.
 *
 * Component layout of the FESystem: components [0, dim) hold the solid
 * velocity and component dim holds the solid volume fraction.
 *
 * @tparam dim Number of spatial dimensions.
 */
template <int dim>
class ContinuumSolidDynamics : public AuxiliaryPhysics<dim, GlobalVectorType>
{
public:
  /**
   * @brief Constructor of the solid phase physics.
   *
   * @param[in] multiphysics_interface Multiphysics interface through which the
   * solid phase exchanges fields with the other physics.
   * @param[in] p_simulation_parameters Simulation parameters.
   * @param[in] p_triangulation Distributed triangulation shared by all physics.
   * @param[in] p_simulation_control Simulation control (time stepping).
   */
  ContinuumSolidDynamics(
    MultiphysicsInterface<dim>      *multiphysics_interface,
    const SimulationParameters<dim> &p_simulation_parameters,
    const std::shared_ptr<parallel::DistributedTriangulationBase<dim>>
                                             &p_triangulation,
    const std::shared_ptr<SimulationControl> &p_simulation_control);

  /**
   * @brief Gather the solid velocity and solid volume fraction for output.
   *
   * @return Vector of output structures holding the solid phase solution.
   */
  std::vector<OutputStruct<dim, GlobalVectorType>>
  gather_output_hook() override;

  /**
   * @brief Finish the simulation (print timers).
   */
  void
  finish_simulation() override;

  /**
   * @brief Shift the solution history for the BDF time-stepping schemes.
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Post-process the solid phase solution. Nothing is computed yet.
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
   * yet for the solid phase.
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
   * @brief Distribute the DoFs, allocate the vectors and matrix, and register
   * the solid phase fields in the multiphysics interface.
   */
  void
  setup_dofs() override;

  /**
   * @brief Set the initial conditions. Fields are set to zero for now.
   */
  void
  set_initial_conditions() override;

  /**
   * @brief Solve the linear system. Not implemented yet.
   */
  void
  solve_linear_system() override;

  /**
   * @brief Indicate whether the solid phase should be solved. Returns false
   * until the governing equations are implemented.
   *
   * @return false.
   */
  bool
  should_solve_auxiliary_physics() override
  {
    return false;
  }

  /**
   * @brief Get the DoFHandler of the solid phase.
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
   * @brief Get the Newton update.
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
   * @brief Get the system right-hand side.
   *
   * @return Reference to the right-hand side vector.
   */
  GlobalVectorType &
  get_system_rhs() override
  {
    return system_rhs;
  }

  /**
   * @brief Get the non-zero constraints.
   *
   * @return Reference to the non-zero constraints.
   */
  AffineConstraints<double> &
  get_nonzero_constraints() override
  {
    return nonzero_constraints;
  }

  /**
   * @brief Print the norms of the Newton update.
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
   * @brief Assemble the system matrix. Not implemented yet.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the system right-hand side. Not implemented yet.
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Define the constraints (hanging nodes only for now).
   */
  void
  define_constraints();

  /// Multiphysics interface used to exchange fields with other physics
  MultiphysicsInterface<dim> *multiphysics;

  /// Timer of the solid phase physics
  TimerOutput computing_timer;

  /// Simulation parameters
  const SimulationParameters<dim> &simulation_parameters;

  /// Triangulation shared by all physics
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation;

  /// Simulation control
  std::shared_ptr<SimulationControl> simulation_control;

  /// DoFHandler of the solid phase
  std::shared_ptr<DoFHandler<dim>> dof_handler;

  /// FESystem: dim velocity components + 1 volume fraction component
  std::shared_ptr<FiniteElement<dim>> fe;

  /// Mapping, cell and face quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::shared_ptr<Quadrature<dim>>     cell_quadrature;
  std::shared_ptr<Quadrature<dim - 1>> face_quadrature;

  /// Extractors of the solid velocity and solid volume fraction
  const FEValuesExtractors::Vector solid_velocity_extractor;
  const FEValuesExtractors::Scalar solid_volume_fraction_extractor;

  /// Locally owned and locally relevant DoFs
  IndexSet locally_owned_dofs;
  IndexSet locally_relevant_dofs;

  /// Solution and solver vectors
  GlobalVectorType                  evaluation_point;
  GlobalVectorType                  local_evaluation_point;
  GlobalVectorType                  newton_update;
  std::shared_ptr<GlobalVectorType> present_solution;
  GlobalVectorType                  system_rhs;
  AffineConstraints<double>         nonzero_constraints;
  AffineConstraints<double>         zero_constraints;
  TrilinosWrappers::SparseMatrix    system_matrix;

  /// Previous solutions for the BDF time-stepping schemes
  std::shared_ptr<std::vector<GlobalVectorType>> previous_solutions;

  /// Solution transfer for mesh adaptation and checkpointing
  std::shared_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;
  std::vector<SolutionTransfer<dim, GlobalVectorType>>
    previous_solutions_transfer;
};

#endif
