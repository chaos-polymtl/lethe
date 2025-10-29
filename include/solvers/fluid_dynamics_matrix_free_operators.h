// SPDX-FileCopyrightText: Copyright (c) 2023-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_matrix_free_operators_h
#define lethe_fluid_dynamics_matrix_free_operators_h

#include <core/mortar_coupling_manager.h>
#include <core/simulation_control.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/base/timer.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

using namespace dealii;

/**
 * @brief Evaluate the value of a function at a batch of points to obtain a vectorized array of numbers
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam Number Abstract type for number across the class (i.e., double).
 * @param function Function to evaluate.
 * @param p_vectorized Batch of points to evaluate function at.
 * @return VectorizedArray<Number> Batch of evaluated values.
 */
template <int dim, typename Number>
VectorizedArray<Number>
evaluate_function(const Function<dim>                       &function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  VectorizedArray<Number> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      result[v] = function.value(p);
    }
  return result;
}

/**
 * @brief Evaluate the gradient of a function at a batch of points to obtain a tensor of vectorized arrays
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam Number Abstract type for number across the class (i.e., double).
 * @param function Function to evaluate.
 * @param p_vectorized Batch of points to evaluate function at.
 * @return Tensor<1, components, VectorizedArray<Number>> Batch of evaluated gradients.
 */
template <int dim, typename Number>
Tensor<1, dim, VectorizedArray<Number>>
evaluate_function_gradient(
  const Function<dim>                       &function,
  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  Tensor<1, dim, VectorizedArray<Number>> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];

      Tensor<1, dim> gradient = function.gradient(p);
      for (int d = 0; d < dim; ++d)
        result[d][v] = gradient[d];
    }
  return result;
}

/**
 * @brief Evaluate the value of a function at a batch of points to obtain a tensor of vectorized arrays
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam Number Abstract type for number across the class (i.e., double).
 * @tparam components Number of solution components.
 * @param function Function to evaluate.
 * @param p_vectorized Batch of points to evaluate function at.
 * @return Tensor<1, components, VectorizedArray<Number>> Batch of evaluated values.
 */
template <int dim, typename Number, int components>
Tensor<1, components, VectorizedArray<Number>>
evaluate_function(const Function<dim>                       &function,
                  const Point<dim, VectorizedArray<Number>> &p_vectorized)
{
  Tensor<1, components, VectorizedArray<Number>> result;
  for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v)
    {
      Point<dim> p;
      for (int d = 0; d < dim; ++d)
        p[d] = p_vectorized[d][v];
      for (unsigned int d = 0; d < components; ++d)
        result[d][v] = function.value(p, d);
    }
  return result;
}

/**
 * @brief A class that serves as base for all the matrix-free
 * Navier-Stokes operators.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesOperatorBase : public EnableObserverPointer
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;
  using value_type       = number;
  using size_type        = VectorizedArray<number>;
  using StabilizationType =
    Parameters::Stabilization::NavierStokesStabilization;

  /**
   * @brief Default constructor.
   */
  NavierStokesOperatorBase();

  /**
   * @brief Constructor that points to the reinit function.
   *
   * @param[in] mapping  Describes the transformations from unit to real cell.
   * @param[in] dof_handler Describes the layout of DoFs and the type of FE.
   * @param[in] constraints Object with constraints according to DoFs.
   * @param[in] quadrature Required for local operations on cells.
   * @param[in] forcing_function Function specified in parameter file as source
   * term.
   * @param[in] properties_manager The physical properties manager (see
   physical_properties_manager.h)
   * @param[in] stabilization Stabilization type specified in parameter file.
   * @param[in] mg_level Level of the operator in case of MG methods.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] physical_properties_manager Required to have the updated values
   * of physical properties in case of ramp or viscous initial conditions.
   * @param[in] boundary_conditions Contains information regarding all boundary
   * conditions. Required to weakly impose boundary conditions.
   * @param[in] enable_hessians_jacobian Flag to turn hessian terms from
   * jacobian on or off.
   * @param[in] enable_hessians_residual Flag to turn hessian terms from
   * residual on or off.
   * @param[in] enable_mortar Flag to turn mortar parameters on or off.

   */
  NavierStokesOperatorBase(
    const Mapping<dim>                  &mapping,
    const DoFHandler<dim>               &dof_handler,
    const AffineConstraints<number>     &constraints,
    const Quadrature<dim>               &quadrature,
    const std::shared_ptr<Function<dim>> forcing_function,
    const std::shared_ptr<PhysicalPropertiesManager>
                                             &physical_properties_manager,
    const StabilizationType                   stabilization,
    const unsigned int                        mg_level,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
    const bool &enable_hessians_jacobian,
    const bool &enable_hessians_residual,
    const bool &enable_mortar);

  /**
   * @brief Initialize the main matrix free object that contains all data and is
   * needed to perform loops over cells, and initialize relevant member
   * variables such as parameters and element size. Calculate edge constraints
   * in case the local smoothing multigrid is used.
   *
   * @param[in] mapping Describes the transformations from unit to real cell.
   * @param[in] dof_handler Describes the layout of DoFs and the type of FE.
   * @param[in] constraints Object with constraints according to DoFs.
   * @param[in] quadrature Required for local operations on cells.
   * @param[in] forcing_function Function specified in parameter file as source
   * term.
   * @param[in] properties_manager The physical properties manager (see
   * physical_properties_manager.h)
   * @param[in] stabilization Stabilization type specified in parameter file.
   * @param[in] mg_level Level of the operator in case of MG methods.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] boundary_conditions Contains information regarding all
   * boundary conditions. Required to weakly impose boundary conditions.
   * @param[in] enable_hessians_jacobian Flag to turn hessian terms from
   * jacobian on or off.
   * @param[in] enable_hessians_residual Flag to turn hessian terms from
   * residual on or off.
   * @param[in] enable_mortar Flag to turn mortar parameters on or off.
   */
  void
  reinit(
    const Mapping<dim>                  &mapping,
    const DoFHandler<dim>               &dof_handler,
    const AffineConstraints<number>     &constraints,
    const Quadrature<dim>               &quadrature,
    const std::shared_ptr<Function<dim>> forcing_function,
    const std::shared_ptr<PhysicalPropertiesManager>
                                             &physical_properties_manager,
    const StabilizationType                   stabilization,
    const unsigned int                        mg_level,
    const std::shared_ptr<SimulationControl> &simulation_control,
    const BoundaryConditions::NSBoundaryConditions<dim> &boundary_conditions,
    const bool &enable_hessians_jacobian,
    const bool &enable_hessians_residual,
    const bool &enable_mortar);

  /**
   * @brief Compute the element size h of the cells required to calculate
   * stabilization parameters.
   */
  void
  compute_element_size();

  /**
   * @brief Precompute forcing term.
   */
  void
  compute_forcing_term();


  /**
   * @brief Precompute buoyancy term for heat-transfer coupling.
   *
   * @param[in] temperature_solution Present solution of the temperature
   * as given by the multiphysics interface.
   * @param[in] temperature_dof_handler DoF Handler used for the heat transfer.
   */
  void
  compute_buoyancy_term(const VectorType      &temperature_solution,
                        const DoFHandler<dim> &temperature_dof_handler);

  /**
   * @brief Get the total number of DoFs.
   *
   * @return Total number of degrees of freedom.
   */
  types::global_dof_index
  m() const;

  /**
   * @brief Access a particular element in the matrix. Only required
   * for compilation and it is not used.
   *
   * @param int
   * @param int
   * @return number
   */
  number
  el(unsigned int, unsigned int) const;

  /**
   * @brief Clear the matrix-free object.
   */
  void
  clear();

  /**
   * @brief Initialize a given vector by delegating it to the MatrixFree
   * function in charge of this task.
   *
   * @param[in,out] vec Vector to be initialized.
   */
  void
  initialize_dof_vector(VectorType &vec) const;

  /**
   * @brief Get the vector partitioner object required for local smoothing
   *  multigrid preconditioner.
   *
   * @return Pointer to vector partitioner.
   */
  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner() const;

  /**
   * @brief Perform an operator evaluation dst = A*src by looping with the help of the
   * MatrixFree object over all cells and evaluating the effect of cell
   * integrals.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Perform the transposed operator evaluation.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Vmult operator for an interface. Required only if local smoothing
   * multigrid is used and for meshes with hanging nodes.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult_interface_down(VectorType &dst, VectorType const &src) const;

  /**
   * @brief Vmult operator for an interface. Required only if local smoothing
   * multigrid is used and for meshes with hanging nodes.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input source vector.
   */
  void
  vmult_interface_up(VectorType &dst, VectorType const &src) const;

  /**
   * @brief Return sparsity pattern used to set up the sparse matrix by
   * get_system_matrix().
   *
   * @return Sparsity pattern.
   */
  const DynamicSparsityPattern &
  get_sparsity_pattern() const;

  /**
   * @brief Calculate matrix if needed, e.g., by coarse-grid solver when a multigrid
   * algorithm is used.
   *
   * @return Trilinos sparse matrix.
   */
  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const;

  /**
   * @brief Get the system matrix free object.
   *
   * @return Matrix free object.
   */
  const MatrixFree<dim, number> &
  get_system_matrix_free() const;

  /**
   * @brief Get the element size object.
   *
   * @return Aligned vector with the size of all elements.
   */
  const AlignedVector<VectorizedArray<number>>
  get_element_size() const;

  /**
   * @brief Compute the diagonal of the operator using an optimized MatrixFree
   * function. Needed for preconditioners.
   *
   * @param[in,out] diagonal The vector where the computed diagonal is stored.
   */
  void
  compute_inverse_diagonal(VectorType &diagonal) const;

  /**
   * @brief Store relevant values of the vector of the last newton step to use it
   * in the Jacobian and pre-calculate the stabilization parameter tau.
   *
   * @param[in] newton_step Vector of the last newton step.
   */
  virtual void
  evaluate_non_linear_term_and_calculate_tau(const VectorType &newton_step);

  /**
   * @brief Store the values of the vector containing the time derivatives of
   * previous solutions to use them in the Jacobian and residual cell integrals.
   *
   * @param[in] time_derivative_previous_solutions Vector with the time
   * derivative of previous solutions.
   */
  void
  evaluate_time_derivative_previous_solutions(
    const VectorType &time_derivative_previous_solutions);

  void
  evaluate_velocity_ale(
    const Mapping<dim>                             &mapping,
    const double                                    radius,
    const Point<dim>                                center_of_rotation,
    std::shared_ptr<Functions::ParsedFunction<dim>> rotor_angular_velocity);

  /**
   * @brief Store the values of the source term calculated if dynamic control
   * is enabled in the appropriate structure.
   *
   * @param[in] beta_force Source term calculated by the dynamic flow control.
   */
  void
  update_beta_force(const Tensor<1, dim> &beta_force);

  /**
   * @brief Evaluate right hand side using the matrix-free operator.
   *
   * @param[in,out] dst Destination vector holding the result.
   * @param[in] src Input vector for which the residual is evaluated.
   */
  void
  evaluate_residual(VectorType &dst, const VectorType &src);

  /**
   * @brief Mortar coupling manager, operator, and evaluator used in the matrix-free
   * solver. The matrix-based solver uses the shared pointers created in the NS
   * base. Because we have a different system_operator here, we need to create
   * new pointers.
   */
  std::shared_ptr<MortarManagerCircle<dim>>      mortar_manager_mf;
  std::shared_ptr<CouplingOperator<dim, double>> mortar_coupling_operator_mf;
  std::shared_ptr<NavierStokesCouplingEvaluation<dim, double>>
    mortar_coupling_evaluator_mf;

protected:
  /**
   * @brief Interface to function that performs a cell integral in a cell batch
   * without gathering and scattering the values. Should be overriden by derived
   * classes.
   *
   * @param[in] integrator FEEvaluation object that allows to evaluate functions
   * at quadrature points and perform cell integrations.
   */
  virtual void
  do_cell_integral_local(FECellIntegrator &integrator) const = 0;

  /**
   * @brief Loop over all cell batches within certain range and perform a cell
   * integral with access to global vectors, i.e., gathering and scattering
   * values.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the cell batch.
   */
  void
  do_cell_integral_range(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const;

  /**
   * @brief Carry out the integration of boundary face integrals.
   *
   * @tparam assemble_residual Flag to assemble the residual or the Jacobian.
   *
   * @param[in] integrator FEFaceEvaluation object that allows to evaluate
   * functions at quadrature points and perform face integrations.
   */
  template <bool assemble_residual>
  void
  do_boundary_face_integral_local(FEFaceIntegrator &integrator) const;

  /**
   * @brief Loop over all boundary face batches within certain range and perform a face
   * integral with access to global vectors, i.e., gathering and scattering
   * values. This is used to weakly impose Dirichlet boundary conditions or
   * outlet boundary conditions.
   *
   * @tparam assemble_residual Flag to assemble the residual or the Jacobian.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the face batch.
   */
  template <bool assemble_residual>
  void
  do_boundary_face_integral_range(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const;

  /**
   * @brief Loop over all internal face batches within certain range and perform a face
   * integral with access to global vectors, i.e., gathering and scattering
   * values. This is only required for compilation and not needed for our CG
   * implementation.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the face batch.
   */
  void
  do_internal_face_integral_range(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const;

  /**
   * @brief Interface to function in charge of computing the residual using a cell
   * integral in a cell batch with gathering and scattering values. Should be
   * overriden by derived classes.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the cell batch.
   */
  virtual void
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const = 0;


private:
  /**
   * @brief Get the refinement edges. Only needed if the local smoothing
   * multigrid approach is used.
   *
   * @param[in] matrix_free Object that contains all data.
   *
   * @return Set containing the indices of the refinement edges.
   */
  static IndexSet
  get_refinement_edges(const MatrixFree<dim, number> &matrix_free);

protected:
  /**
   * @brief Object storing all the relevant data needed for the matrix free
   * implementation.
   *
   */
  MatrixFree<dim, number> matrix_free;

  /**
   * @brief Object that stores constraints for the matrix free implementation.
   *
   */
  AffineConstraints<number> constraints;

  /**
   * @brief Sparsity pattern used when the matrix of certain level is computed
   * and stored.
   */
  mutable DynamicSparsityPattern dsp;

  /**
   * @brief Sparse trilinos matrix used when the matrix of certain level is computed
   * and stored.
   *
   */
  mutable TrilinosWrappers::SparseMatrix system_matrix;

  /**
   * @brief Aligned vector to store the element size of all the elements.
   *
   */
  AlignedVector<VectorizedArray<number>> element_size;

  /**
   * @brief Finite element degree for the operator.
   *
   */
  unsigned int fe_degree;

  /**
   * @brief Force function or source function for the Navier-Stokes equations.
   *
   */
  std::shared_ptr<Function<dim>> forcing_function;

  /**
   * @brief Additional source term in the case of dynamic flow control.
   *
   */
  Tensor<1, dim, VectorizedArray<number>> beta_force;

  /**
   * @brief Stabilization type needed to add or remove terms from operator.
   *
   */
  StabilizationType stabilization;

  /**
   * @brief Object storing the information regarding the time stepping method.
   *
   */
  std::shared_ptr<SimulationControl> simulation_control;

  /**
   * @brief Object pointing to the physical properties manager of the matrix free class.
   *
   */
  std::shared_ptr<PhysicalPropertiesManager> properties_manager;

  /**
   * @brief  Boundary conditions object to impose the correct boundary conditions. This object
   * is only used to impose boundary conditions that appear in the weak form as
   * face terms (Weak Dirichlet boundary conditions and, eventually, outlets).
   * However, the entire boundary_conditions object is stored instead of trying
   * to isolate which parameters are used since many of them are used in the
   * initialization of the boundary conditions.
   *
   */
  BoundaryConditions::NSBoundaryConditions<dim> boundary_conditions;

  /**
   * @brief Flag to turn hessian terms from jacobian on or off.
   *
   */
  bool enable_hessians_jacobian;

  /**
   * @brief Flag to turn hessian terms from residual on or off.
   *
   */
  bool enable_hessians_residual;

  /**
   * @brief Table with precomputed forcing term values.
   *
   */
  Table<2, Tensor<1, dim, VectorizedArray<number>>> forcing_terms;


  /**
   * @brief Table with precomputed buoyancy term values.
   *
   */
  Table<2, Tensor<1, dim, VectorizedArray<number>>> buoyancy_term;

  /**
   * @brief Flag to turn the calculation of face terms on or off.
   * This is used for weakly imposed Dirichlet boundary conditions or outlets.
   *
   */
  bool enable_face_terms;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the previous newton step.
   *
   */
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    nonlinear_previous_values;

  /**
   * @brief Table with correct alignment for vectorization to store the gradients
   * of the previous newton step.
   *
   */
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_gradient;

  /**
   * @brief Table with correct alignment for vectorization to store the hessian
   * diagonal of the previous newton step.
   *
   */
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_hessian_diagonal;

  /**
   * @brief Table with correct alignment for vectorization to store the time derivatives
   * of previous solutions.
   *
   */
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    time_derivatives_previous_solutions;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the velocity at faces.
   *
   */
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    face_nonlinear_previous_values;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the target velocity of a weakly imposed Dirichlet boundary condition.
   *
   */
  Table<2, Tensor<1, dim, VectorizedArray<number>>> face_target_velocity;


  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the stabilization parameter tau.
   *
   */
  Table<2, VectorizedArray<number>> stabilization_parameter;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the stabilization parameter tau lsic.
   *
   */
  Table<2, VectorizedArray<number>> stabilization_parameter_lsic;


  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the kinematic viscosity.
   *
   */
  Table<2, VectorizedArray<number>> kinematic_viscosity_vector;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the gradient of the kinematic viscosity w.r.t shear rate.
   *
   */
  Table<2, VectorizedArray<number>> grad_kinematic_viscosity_shear_rate;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the kinematic viscosity gradient.
   *
   */
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    kinematic_viscosity_gradient;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the shear_rate.
   *
   */
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    previous_shear_rate;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the shear_rate_magnitude.
   *
   */
  Table<2, VectorizedArray<number>> previous_shear_rate_magnitude;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the face penalization term effective_beta_face
   *
   */
  Table<1, VectorizedArray<number>> effective_beta_face;

  /**
   * @brief Flag to turn the computation of mortar coupling terms on or off.
   */
  bool enable_mortar;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the ALE velocity used in mortar coupling terms.
   *
   */
  Table<2, Tensor<1, dim, VectorizedArray<number>>> velocity_ale;

  /**
   * @brief Vector with the constrained indices used for the local smoothing approach.
   *
   */
  std::vector<unsigned int> constrained_indices;

  /**
   * @brief Vector with the constrained values used for the local smoothing approach.
   *
   */
  mutable std::vector<std::pair<number, number>> constrained_values;

  /**
   * @brief Vector with the constrained edge indices for the local smoothing approach.
   *
   */
  std::vector<unsigned int> edge_constrained_indices;

  /**
   * @brief Variable to identify whether there are edge constrained indices for the
   * local smoothing approach.
   *
   */
  bool has_edge_constrained_indices = false;

  /**
   * @brief Vector with edge constrained values used for the local smoothing approach.
   *
   */
  mutable std::vector<number> edge_constrained_values;


  /**
   * @brief Vector used to mark edge constrained cells in order to verify constrained indices
   * for the local smoothing approach.
   *
   */
  std::vector<bool> edge_constrained_cell;

  /**
   * @brief DoF mask object needed for the sparsity pattern and computation of the system
   * matrix of the coarse level in case FE_Q_iso_Q1 elements are used.
   *
   */
  Table<2, bool> bool_dof_mask;

  /**
   * @brief Conditional OStream for parallel output.
   *
   */
  ConditionalOStream pcout;

public:
  /**
   * @brief Timer for internal operator calls.
   *
   */
  mutable TimerOutput timer;
};

/**
 * @brief Implements the matrix-free operator to solve steady/transient Navier-Stokes equations
 * using stabilization.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesStabilizedOperator
  : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;

  NavierStokesStabilizedOperator();

protected:
  /**
   * @brief Perform cell integral on a cell batch without gathering and scattering
   * the values, and according to the Jacobian of the discretized
   * steady/transient Navier-Stokes equations with stabilization.
   *
   * @param[in] integrator FEEvaluation object that allows to evaluate functions
   * at quadrature points and perform cell integrations.
   */
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;

  /**
   * @brief Perform cell integral on a cell batch with gathering and scattering
   * the values, and according to the residual of the discretized
   * steady/transient Navier-Stokes equations with stabilization.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the cell batch.
   */
  void
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const override;
};

/**
 * @brief Implements the matrix-free operator to solve steady/transient non-Newtonian Navier-Stokes equations using stabilization.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesNonNewtonianStabilizedOperator
  : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;

  NavierStokesNonNewtonianStabilizedOperator();

protected:
  /**
   * @brief Perform cell integral on a cell batch without gathering and scattering
   * the values, and according to the Jacobian of the discretized
   * steady/transient Navier-Stokes equations with stabilization.
   *
   * @param[in] integrator FEEvaluation object that allows to evaluate functions
   * at quadrature points and perform cell integrations.
   */
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;

  /**
   * @brief Perform cell integral on a cell batch with gathering and scattering
   * the values, and according to the residual of the discretized
   * steady/transient Navier-Stokes equations with stabilization.
   *
   * @param[in] matrix_free Object that contains all data.
   * @param[in,out] dst Global vector where the final result is added.
   * @param[in] src Input vector with all values in all cells.
   * @param[in] range Range of the cell batch.
   */
  void
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const override;
};

#endif
