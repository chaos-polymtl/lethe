// SPDX-FileCopyrightText: Copyright (c) 2025-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#ifndef lethe_fluid_dynamics_vans_matrix_free_operators_h
#define lethe_fluid_dynamics_vans_matrix_free_operators_h

#include <solvers/fluid_dynamics_matrix_free_operators.h>

#include <fem-dem/particle_projector.h>

using namespace dealii;

/**
 * @brief Implements the matrix-free operator for the VANS equations
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


  /**
   * @brief Constructor for the VANS operator
   *
   * @param[in] cfd_dem_parameters The CFD-DEM parameters used for the
   * simulation. These parameters are notably used to distinguish between
   * model A and model B as well as to establish grad-div stabilization.
   */
  VANSOperator(const Parameters::CFDDEM &cfd_dem_parameters)
    : NavierStokesOperatorBase<dim, number>()
    , cfd_dem_parameters(cfd_dem_parameters)
  {
    AssertThrow(
      this->cfd_dem_parameters.vans_model == Parameters::VANSModel::modelA,
      ExcMessage(
        "Only model A of the VANS equation is supported by the VANS matrix-free solver"));
  }

  /**
   * @brief Compute the void fraction at the fluid quadrature point.
   * This essentially calculates the void fraction and the void fraction
   * gradient.
   *
   * @param[in] void_fraction_dof_handler The dof handler associated with the
   * void fraction
   *
   * @param[in] void_fraction_solution The solution of the void fraction
   *
   */
  void
  compute_void_fraction(
    const DoFHandler<dim>                            &void_fraction_dof_handler,
    const LinearAlgebra::distributed::Vector<double> &void_fraction_solution);

  /**
   * @brief Compute the terms required to assemble the particle-fluid forces at
   * the gauss points. The function gathers the particle-fluid force (without
   * drag), the particle-fluid drag and the particle velocity. With these three
   * elements, the full particle-fluid coupling and its jacobian can be
   * established within the matrix-free operator.
   *
   * @param[in] pf_force_dof_handler The dof handler associated with the
   * particle-fluid forces.
   * @param[in] pf_force_solution The solution of the particle-fluid forces.
   * @param[in] pf_drag_dof_handler The dof handler associated with the
   * drag force.
   * @param[in] pf_drag_solution The solution of the drag force.
   * @param[in] particle_velocity_dof_handler The dof handler associated with
   * the particle velocity.
   * @param[in] pf_drag_solution The solution of the particle velocity.
   *
   */
  void
  compute_particle_fluid_force(
    const DoFHandler<dim>                            &pf_force_dof_handler,
    const LinearAlgebra::distributed::Vector<double> &pf_force_solution,
    const DoFHandler<dim>                            &pf_drag_dof_handler,
    const LinearAlgebra::distributed::Vector<double> &pf_drag_solution,
    const DoFHandler<dim> &particle_velocity_dof_handler,
    const LinearAlgebra::distributed::Vector<double>
      &particle_velocity_solution);

protected:
  /**
   * @brief Store relevant values of the vector of the last newton step to use it
   * in the Jacobian and pre-calculate the stabilization parameters tau and
   * gamma.
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

  /// Table with correct alignment for vectorization to store the values of the
  /// void fraction
  Table<2, VectorizedArray<number>> void_fraction;

  /// Table with correct alignment for vectorization to store the values of the
  /// void fraction gradient
  Table<2, Tensor<1, dim, VectorizedArray<number>>> void_fraction_gradient;

  /// Table with correct alignment for vectorization to store the values of the
  /// particle-fluid forces (without drag). These are the resulting forces
  /// applied on the particles due to hydrodynamic interaction, normalized by
  /// the volume of the cells.
  Table<2, Tensor<1, dim, VectorizedArray<number>>> particle_fluid_force;

  /// Table with correct alignment for vectorization to store the values of the
  /// particle-fluid drag.
  Table<2, Tensor<1, dim, VectorizedArray<number>>> particle_fluid_drag;

  /// Table with correct alignment for vectorization to store the values of the
  /// particle-fluid drag.
  Table<2, Tensor<1, dim, VectorizedArray<number>>> particle_velocity;

  /// Table with correct alignment for vectorization to store the values of the
  /// grad-div gamma parameter
  Table<2, VectorizedArray<number>> grad_div_gamma;

  ///  Internal copy of the CFD-DEM parameters. This is used for grad-div
  ///  stabilization, but also to switch between form A and B of the VANS
  ///  equations.
  const Parameters::CFDDEM cfd_dem_parameters;
};

#endif
