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

#include <solvers/simulation_parameters.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

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

  using VectorType = LinearAlgebra::distributed::Vector<number>;

  using value_type = number;

  using size_type = VectorizedArray<number>;

  /**
   * @brief Construct a new Navier Stokes Operator Base object.
   */
  NavierStokesOperatorBase();

  /**
   * @brief Construct a new Navier Stokes Operator Base object and point to
   * the reinit function of this class.
   *
   * @param mapping  Describes the transformations from unit to real cell.
   * @param dof_handler Describes the layout of DoFs and the type of FE.
   * @param constraints Object with constraints according to DoFs.
   * @param quadrature Required for local operations on cells.
   * @param forcing_function Function specified in parameter file as source term.
   * @param kinematic_viscosity Kinematic viscosity.
   * @param mg_level Level of the operator in case of MG methods.
   */
  NavierStokesOperatorBase(
    const Mapping<dim>                &mapping,
    const DoFHandler<dim>             &dof_handler,
    const AffineConstraints<number>   &constraints,
    const Quadrature<dim>             &quadrature,
    const Function<dim>               *forcing_function,
    const double                       kinematic_viscosity,
    const unsigned int                 mg_level,
    std::shared_ptr<SimulationControl> simulation_control);
  /**
   * @brief Initialize the main matrix free object that contains all data and is
   * needed to perform loops over cells, and initialize relevant member
   * variables such as parameters and element size.
   *
   * @param mapping Describes the transformations from unit to real cell.
   * @param dof_handler Describes the layout of DoFs and the type of FE.
   * @param constraints Object with constraints according to DoFs.
   * @param quadrature Required for local operations on cells.
   * @param forcing_function Function specified in parameter file as source term.
   * @param kinematic_viscosity Kinematic viscosity.
   * @param mg_level Level of the operator in case of MG methods.
   */
  void
  reinit(const Mapping<dim>                &mapping,
         const DoFHandler<dim>             &dof_handler,
         const AffineConstraints<number>   &constraints,
         const Quadrature<dim>             &quadrature,
         const Function<dim>               *forcing_function,
         const double                       kinematic_viscosity,
         const unsigned int                 mg_level,
         std::shared_ptr<SimulationControl> simulation_control);

  /**
   * @brief Compute the element size h of the cells required to calculate
   * stabilization parameters.
   */
  void
  compute_element_size();

  /**
   * @brief Return the total number of DoFs.
   *
   * @return types::global_dof_index.
   */
  types::global_dof_index
  m() const;

  /**
   * @brief Access a particular element in the matrix. This function is only required
   * for compilation and it is not used.
   *
   * @param int
   * @param int
   * @return number
   */
  number
  el(unsigned int, unsigned int) const;

  /**
   * @brief This function is only required for compilation and it is not used.
   */
  void
  clear();

  /**
   * @brief Initialize a given vector by delegating it to the MatrixFree
   * function in charge of this task.
   *
   * @param vec Vector to be initialized.
   */
  void
  initialize_dof_vector(VectorType &vec) const;

  /**
   * @brief Get the vector partitioner object required for LS multigrid
   */
  const std::shared_ptr<const Utilities::MPI::Partitioner> &
  get_vector_partitioner() const;

  /**
   * @brief Perform an operator evaluation dst = A*src by looping with the help of the
   * MatrixFree object over all cells and evaluatinf the effect of cell
   * integrals.
   *
   * @param dst Destination vector holding the result.
   * @param src Input vector.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Perform the transposed operator evaluation
   *
   * @param dst Destination vector holding the result.
   * @param src Input vector.
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Vmult ooperator for an interface. Required only if local smoothing
   * multigrid is used and for meshes with hanging nodes.
   *
   * @param dst Destination vector holding the result.
   * @param src Input vector.
   */
  void
  vmult_interface_down(VectorType &dst, VectorType const &src) const;

  /**
   * @brief Vmult ooperator for an interface. Required only if local smoothing
   * multigrid is used and for meshes with hanging nodes.
   *
   * @param dst Destination vector holding the result.
   * @param src Input vector.
   */
  void
  vmult_interface_up(VectorType &dst, VectorType const &src) const;

  /**
   * @brief Calculate matrix if needed, e.g., by coarse-grid solver when a multigrid
   * algorithm is used.
   *
   * @return const TrilinosWrappers::SparseMatrix&.
   */
  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const;

  /**
   * @brief Get the system matrix free object.
   *
   * @return const MatrixFree<dim, number>& .
   */
  const MatrixFree<dim, number> &
  get_system_matrix_free() const;

  /**
   * @brief Get the element size object.
   *
   * @return const AlignedVector<VectorizedArray<number>>.
   */
  const AlignedVector<VectorizedArray<number>>
  get_element_size() const;

  /**
   * @brief Compute the diagonal of the operator using an optimized MatrixFree
   * function. Needed for preconditioners.
   *
   * @param diagonal The vector where the computed diagonal is stored.
   */
  void
  compute_inverse_diagonal(VectorType &diagonal) const;

  /**
   * @brief Store relevant values of the vector of the last newton step to use it
   * in the Jacobian.
   *
   * @param newton_step Vector of the last newton step.
   */
  void
  evaluate_non_linear_term(const VectorType &newton_step);

  /**
   * @brief Store the values of the vector containing the time derivatives of
   * previous solutions to use them in the Jacobian and residual cell integrals
   *
   * @param time_derivative_previous_solutions Vector with the time derivative
   * of previous solutions.
   */
  void
  evaluate_time_derivative_previous_solutions(
    const VectorType &time_derivative_previous_solutions);

  /**
   * @brief Evaluate right hand side using the matrix-free operator
   *
   * @param dst Destination vector holding the result
   * @param src Input vector for which the residual is evaluated
   */
  void
  evaluate_residual(VectorType &dst, const VectorType &src);

  /**
   * @brief Sets the kinematic viscosity in the operator
   *
   * @param p_kinematic_viscosity New value of the kinematic viscosity
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
   * @param integrator FEEvaluation object that allows to evaluate functions at
   * quadrature points and perform cell integrations.
   */
  virtual void
  do_cell_integral_local(FECellIntegrator &integrator) const = 0;

  /**
   * @brief Loop over all cell batches withing certain range and perform a cell
   * integral with access to global vectors, i.e., gathering and scattering
   * values.
   *
   * @param matrix_free Object that contains all data.
   * @param dst Global vector where the final result is added.
   * @param src Input vector with all values in all cells.
   * @param range Range of the cell batch.
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
   * @param matrix_free Object that contains all data.
   * @param dst Global vector where the final result is added.
   * @param src Input vector with all values in all cells.
   * @param range Range of the cell batch.
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
   * @param matrix_free Object that contains all data.
   * @return IndexSet Set containing the indices of the refinement edges.
   */
  static IndexSet
  get_refinement_edges(const MatrixFree<dim, number> &matrix_free);

protected:
  MatrixFree<dim, number>                matrix_free;
  AffineConstraints<number>              constraints;
  mutable TrilinosWrappers::SparseMatrix system_matrix;
  AlignedVector<VectorizedArray<number>> element_size;
  unsigned int                           fe_degree;
  const Function<dim>                   *forcing_function;
  double                                 kinematic_viscosity;
  std::shared_ptr<SimulationControl>     simulation_control;

  // Variables needed from the last Newton step vector
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    nonlinear_previous_values;
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_gradient;
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_hessian_diagonal;

  // Variable needed to store the time derivative of the previous solutions
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    time_derivatives_previous_solutions;

  // Variables needed for the local smoothing approach
  std::vector<unsigned int>                      constrained_indices;
  mutable std::vector<std::pair<number, number>> constrained_values;
  std::vector<unsigned int>                      edge_constrained_indices;
  bool has_edge_constrained_indices = false;
  mutable std::vector<std::pair<number, number>> edge_constrained_values;
  std::vector<bool>                              edge_constrained_cell;
};

/**
 * @brief Class in charge of implementing the main function required to
 * solve the Navier-Stokes equations using SUPG/PSPG stabilization and the
 * matrix-free approach.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesSUPGPSPGOperator
  : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;

  NavierStokesSUPGPSPGOperator();

protected:
  /**
   * @brief Perform cell integral on a cell batch without gathering and scattering
   * the values, and according to the Jacobian of the Navier-Stokes equations
   * with SUPG/PSPG stabilization.
   *
   * @param integrator FEEvaluation object that allows to evaluate functions at
   * quadrature points and perform cell integrations.
   */
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;

  /**
   * @brief Perform cell integral on a cell batch with gathering and scattering
   * the values, and according to the residual of the Navier-Stokes equations
   * with SUPG/PSPG stabilization.
   *
   * @param matrix_free Object that contains all data.
   * @param dst Global vector where the final result is added.
   * @param src Input vector with all values in all cells.
   * @param range Range of the cell batch.
   */
  void
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const override;
};

/**
 * @brief Class in charge of implementing the main function required to
 * solve the transient Navier-Stokes equations using SUPG/PSPG stabilization
 * and the matrix-free approach.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class NavierStokesTransientSUPGPSPGOperator
  : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;

  NavierStokesTransientSUPGPSPGOperator();

protected:
  /**
   * @brief Perform cell integral on a cell batch without gathering and scattering
   * the values, and according to the Jacobian of the Navier-Stokes equations
   * with SUPG/PSPG stabilization.
   *
   * @param integrator FEEvaluation object that allows to evaluate functions at
   * quadrature points and perform cell integrations.
   */
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;

  /**
   * @brief Perform cell integral on a cell batch with gathering and scattering
   * the values, and according to the residual of the Navier-Stokes equations
   * with SUPG/PSPG stabilization.
   *
   * @param matrix_free Object that contains all data.
   * @param dst Global vector where the final result is added.
   * @param src Input vector with all values in all cells.
   * @param range Range of the cell batch.
   */
  void
  local_evaluate_residual(
    const MatrixFree<dim, number>               &matrix_free,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::pair<unsigned int, unsigned int> &range) const override;
};


#endif
