// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_matrix_free_operators_h
#define lethe_fluid_dynamics_vans_matrix_free_operators_h

#include <solvers/fluid_dynamics_matrix_free_operators.h>

#include <fem-dem/void_fraction.h>

using namespace dealii;

/**
 * @brief A class that serves as base for all the matrix-free
 * Navier-Stokes operators.
 *
 * @tparam dim An integer that denotes the number of spatial dimensions.
 * @tparam number Abstract type for number across the class (i.e., double).
 */
template <int dim, typename number>
class VANSOperator : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  using FEFaceIntegrator = FEFaceEvaluation<dim, -1, 0, dim + 1, number>;
  using VectorType       = LinearAlgebra::distributed::Vector<number>;


  // TODO - Document
  VANSOperator();

  // TODO - Document
  void
  evaluate_void_fraction(const VoidFractionBase<dim> &void_fraction_manager);

protected:
  /**
   * @brief Store relevant values of the vector of the last newton step to use it
   * in the Jacobian and pre-calculate the stabilization parameter tau.
   *
   * @param[in] newton_step Vector of the last newton step.
   */
  virtual void
  evaluate_non_linear_term_and_calculate_tau(
    const VectorType &newton_step) override;



  /**
   * @brief Perform cell integral on a cell batch without gathering and scattering
   * the values, and according to the Jacobian of the discretized
   * steady/transient VANS equations with stabilization.
   *
   * @param[in] integrator FEEvaluation object that allows to evaluate functions
   * at quadrature points and perform cell integrations.
   */
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;

  /**
   * @brief Perform cell integral on a cell batch with gathering and scattering
   * the values, and according to the residual of the discretized
   * steady/transient VANS equations with stabilization.
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


  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the void fraction
   *
   */
  Table<2, VectorizedArray<number>> void_fraction;

  /**
   * @brief Table with correct alignment for vectorization to store the values
   * of the void fraction gradient
   *
   */
  Table<2, Tensor<1, dim, VectorizedArray<number>>> void_fraction_gradient;
};

#endif
