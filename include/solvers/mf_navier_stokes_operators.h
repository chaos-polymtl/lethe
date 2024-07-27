/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 3.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------*/

#ifndef lethe_mf_navier_stokes_operators_h
#define lethe_mf_navier_stokes_operators_h

#include <core/bdf.h>
#include <core/simulation_control.h>
#include <core/time_integration_utilities.h>

#include <solvers/simulation_parameters.h>

#include <deal.II/base/timer.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

#include <deal.II/multigrid/mg_tools.h>

using namespace dealii;

/**
 * @brief A class that serves as base for all the matrix-free
 * Navier-Stokes operators.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesOperatorBase : public Subscriptor
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
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
   * @param[in] kinematic_viscosity Kinematic viscosity.
   * @param[in] stabilization Stabilization type specified in parameter file.
   * @param[in] mg_level Level of the operator in case of MG methods.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] enable_hessians_jacobian Flag to turn hessian terms from
   * jacobian on or off.
   * @param[in] enable_hessians_residual Flag to turn hessian terms from
   * residual on or off.
   */
  NavierStokesOperatorBase(
    const Mapping<dim>                  &mapping,
    const DoFHandler<dim>               &dof_handler,
    const AffineConstraints<number>     &constraints,
    const Quadrature<dim>               &quadrature,
    const std::shared_ptr<Function<dim>> forcing_function,
    const double                         kinematic_viscosity,
    const StabilizationType              stabilization,
    const unsigned int                   mg_level,
    std::shared_ptr<SimulationControl>   simulation_control,
    const bool                          &enable_hessians_jacobian,
    const bool                          &enable_hessians_residual);

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
   * @param[in] kinematic_viscosity Kinematic viscosity.
   * @param[in] stabilization Stabilization type specified in parameter file.
   * @param[in] mg_level Level of the operator in case of MG methods.
   * @param[in] simulation_control Required to get the time stepping method.
   * @param[in] enable_hessians_jacobian Flag to turn hessian terms from
   * jacobian on or off.
   * @param[in] enable_hessians_residual Flag to turn hessian terms from
   * residual on or off.
   */
  void
  reinit(const Mapping<dim>                  &mapping,
         const DoFHandler<dim>               &dof_handler,
         const AffineConstraints<number>     &constraints,
         const Quadrature<dim>               &quadrature,
         const std::shared_ptr<Function<dim>> forcing_function,
         const double                         kinematic_viscosity,
         const StabilizationType              stabilization,
         const unsigned int                   mg_level,
         std::shared_ptr<SimulationControl>   simulation_control,
         const bool                          &enable_hessians_jacobian,
         const bool                          &enable_hessians_residual);

  /**
   * @brief Compute the element size h of the cells required to calculate
   * stabilization parameters.
   */
  void
  compute_element_size();

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
  void
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
   * @brief Sets the kinematic viscosity in the operator.
   *
   * @param[in] p_kinematic_viscosity New value of the kinematic viscosity.
   */
  void
  set_kinematic_viscosity(const double p_kinematic_viscosity)
  {
    kinematic_viscosity = p_kinematic_viscosity;
  }

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
   * @brief Kinematic viscosity needed for the operator.
   *
   */
  double kinematic_viscosity;


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

#endif
