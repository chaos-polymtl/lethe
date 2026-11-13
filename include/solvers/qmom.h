// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

/**
 * @file qmom.h
 * @brief Auxiliary physics for the transport of a particle size distribution
 * with the quadrature method of moments (QMOM), coupled to the QMOMKit
 * library.
 *
 * The unknowns are regular scalars and the first raw moments of the particle
 * size distribution, called the fields of the QMOM physics. All the fields
 * are discretized with the same scalar discontinuous Galerkin element, on a
 * single DoFHandler, and each field is stored in its own vector. Since all
 * the fields are advected by the same fluid velocity, the matrix of their
 * transport is the same: it is assembled and preconditioned once per time
 * step, then solved with one right-hand side per field.
 *
 * Each time step of the fluid dynamics is split (Lie splitting) into:
 *  - the transport of all the fields in physical space by the fluid velocity,
 *    solved with an upwind discontinuous Galerkin discretization and an
 *    implicit BDF time integration (advance_physical_space());
 *  - the evolution of the moments in phase space (nucleation, growth and
 *    coalescence), applied independently at every degree of freedom
 *    (advance_phase_space()).
 *
 * The transport is linear in the fields, since the fluid velocity is given.
 * The QMOM physics is therefore driven by the linear solver strategy:
 * assemble_system_matrix() assembles the transport matrix,
 * setup_preconditioner() builds its preconditioner, assemble_system_rhs()
 * assembles the right-hand side of every field and solve_linear_system()
 * advances the fields over the time step.
 */

#ifndef lethe_qmom_h
#define lethe_qmom_h

#include <core/boundary_conditions.h>
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
 * Field layout: field i < n_scalars is the scalar i and field n_scalars + k is
 * the raw moment M_k of the particle size distribution, k < n_moments. This is
 * the ordering of the components of the QMOMFieldFunction of the initial
 * conditions, boundary conditions and analytical solution.
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
   * @brief Gather the fields for output, one output structure per field.
   *
   * @return Vector of output structures holding the QMOM fields.
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
   * @brief Shift the solution history of every field.
   */
  void
  percolate_time_vectors() override;

  /**
   * @brief Post-process the QMOM fields: compute the L2 error of every field
   * if the analytical solution is enabled.
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
   * @brief Prepare the present and previous solutions of every field for mesh
   * refinement/coarsening.
   */
  void
  pre_mesh_adaptation() override;

  /**
   * @brief Interpolate the present and previous solutions of every field onto
   * the adapted mesh.
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
   * @brief Serialize the present and previous solutions of every field for
   * restart.
   */
  void
  write_checkpoint() override;

  /**
   * @brief Deserialize the present and previous solutions of every field.
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
   * in the multiphysics interface. The fields themselves are not shared with
   * the other physics.
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
   * @brief Update the time of the Dirichlet boundary conditions if they are
   * time-dependent.
   */
  void
  update_boundary_conditions() override;

  /**
   * @brief Build the ILU preconditioner of the transport matrix, shared by
   * all the fields.
   */
  void
  setup_preconditioner() override;

  /**
   * @brief Advance the fields over the current time step with Lie splitting:
   * transport in physical space, then evolution in phase space. This is the
   * last function called by the linear solver strategy.
   */
  void
  solve_linear_system() override;

  /**
   * @brief Get the DoFHandler shared by all the QMOM fields.
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
   * field.
   *
   * @return Never returns.
   */
  GlobalVectorType &
  get_present_solution() override
  {
    AssertThrow(false, ExcMessage(single_vector_error_message));
  }

  /**
   * @brief Not available: the QMOM physics has one right-hand side per field.
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
   * @brief Assemble the transport matrix shared by all the fields, with the
   * velocity of the fluid dynamics solution.
   */
  void
  assemble_system_matrix() override;

  /**
   * @brief Assemble the right-hand side of every field, with the velocity of
   * the fluid dynamics solution.
   */
  void
  assemble_system_rhs() override;

  /**
   * @brief Assemble the transport matrix: BDF time derivative, upwind
   * discontinuous Galerkin advection and outflow through the boundaries.
   *
   * @tparam FluidVectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_solution Ghosted fluid dynamics solution.
   */
  template <typename FluidVectorType>
  void
  assemble_transport_matrix(const FluidVectorType &fluid_solution);

  /**
   * @brief Assemble the right-hand side of every field: previous solutions of
   * the BDF time derivative and inflow through the Dirichlet boundaries.
   *
   * @tparam FluidVectorType Vector type of the fluid dynamics solution.
   * @param[in] fluid_solution Ghosted fluid dynamics solution.
   */
  template <typename FluidVectorType>
  void
  assemble_transport_rhs(const FluidVectorType &fluid_solution);

  /**
   * @brief Get the type of the QMOM boundary condition of a boundary.
   *
   * @param[in] boundary_id Boundary id.
   *
   * @return Type of the boundary condition, none if the boundary has no QMOM
   * boundary condition.
   */
  BoundaryConditions::BoundaryType
  get_boundary_type(const types::boundary_id boundary_id) const;

  /**
   * @brief Define the (empty) constraints of the discontinuous Galerkin
   * discretization.
   */
  void
  define_constraints();

  /**
   * @brief Transport all the fields in physical space over the time step:
   * solve the transport system of every field, with the matrix and the
   * right-hand sides assembled by the linear solver strategy.
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
   * @brief Compute the L2 error of every field with respect to the analytical
   * solution at the current time.
   *
   * @return L2 errors, ordered as the fields.
   */
  std::vector<double>
  calculate_L2_errors();

  /**
   * @brief Gather pointers to the present solutions of all the fields,
   * followed by their previous solutions, in the order used by the
   * checkpoints and the mesh adaptation.
   *
   * @return Pointers to all the solution vectors.
   */
  std::vector<GlobalVectorType *>
  gather_all_solutions();

  /**
   * @brief Transfer all the solution vectors, ordered as in
   * gather_all_solutions(), after a mesh adaptation or from a checkpoint.
   *
   * @tparam TransferFunction Callable taking a std::vector<GlobalVectorType *>
   * of locally owned vectors and filling them (with
   * SolutionTransfer::interpolate() or SolutionTransfer::deserialize()).
   * @param[in] transfer Function that fills the locally owned vectors.
   */
  template <typename TransferFunction>
  void
  transfer_all_solutions(const TransferFunction &transfer);

  /// Error message of the PhysicsSolver functions that require a single
  /// solution vector
  static constexpr const char *single_vector_error_message =
    "The QMOM physics stores one solution vector per field and is solved with "
    "the linear solver strategy: this function should not be called.";

  /// Number of regular scalars
  const unsigned int n_scalars;

  /// Number of transported moments
  const unsigned int n_moments;

  /// Total number of fields: n_scalars + n_moments
  const unsigned int n_fields;

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

  /// DoFHandler shared by all the fields
  std::shared_ptr<DoFHandler<dim>> dof_handler;

  /// Scalar discontinuous element shared by all the fields
  std::unique_ptr<FiniteElement<dim>> fe;

  /// Mapping, cell and face quadratures
  std::shared_ptr<Mapping<dim>>        mapping;
  std::unique_ptr<Quadrature<dim>>     cell_quadrature;
  std::unique_ptr<Quadrature<dim - 1>> face_quadrature;

  /// Names of the fields
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

  /// Empty constraints of the discontinuous Galerkin discretization
  AffineConstraints<double> constraints;

  /// Transport matrix, shared by all the fields
  TrilinosWrappers::SparseMatrix system_matrix;

  /// ILU preconditioner of the transport matrix
  std::shared_ptr<TrilinosWrappers::PreconditionILU> ilu_preconditioner;

  /// Ghosted present solution of each field: present_solutions[field]
  std::vector<GlobalVectorType> present_solutions;

  /// Ghosted previous solutions of each field, used by the BDF schemes:
  /// previous_solutions[field][p] is the solution p + 1 time steps ago
  std::vector<std::vector<GlobalVectorType>> previous_solutions;

  /// Right-hand side of each field: system_rhs[field]
  std::vector<GlobalVectorType> system_rhs;

  /// Solution transfer of all the solution vectors, for mesh adaptation and
  /// checkpointing
  std::unique_ptr<SolutionTransfer<dim, GlobalVectorType>> solution_transfer;

  /// L2 errors of the fields with respect to the analytical solution
  ConvergenceTable error_table;
};

#endif
