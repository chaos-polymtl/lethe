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

#include <solvers/simulation_parameters.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

using namespace dealii;

template <int dim, typename number>
class NavierStokesOperatorBase : public Subscriptor
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;

  using VectorType = LinearAlgebra::distributed::Vector<number>;

  using value_type = number;

  using size_type = VectorizedArray<number>;

  NavierStokesOperatorBase();

  NavierStokesOperatorBase(const Mapping<dim> &             mapping,
                           const DoFHandler<dim> &          dof_handler,
                           const AffineConstraints<number> &constraints,
                           const Quadrature<dim> &          quadrature,
                           const SimulationParameters<dim> &nsparam,
                           const unsigned int               mg_level);

  void
  reinit(const Mapping<dim> &             mapping,
         const DoFHandler<dim> &          dof_handler,
         const AffineConstraints<number> &constraints,
         const Quadrature<dim> &          quadrature,
         const SimulationParameters<dim> &nsparam,
         const unsigned int               mg_level);

  void
  compute_element_size();

  types::global_dof_index
  m() const;

  number
  el(unsigned int, unsigned int) const;

  void
  clear();

  void
  initialize_dof_vector(VectorType &vec) const;

  void
  vmult(VectorType &dst, const VectorType &src) const;

  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  void
  vmult_interface_down(VectorType &dst, VectorType const &src) const;

  void
  vmult_interface_up(VectorType &dst, VectorType const &src) const;

  const TrilinosWrappers::SparseMatrix &
  get_system_matrix() const;

  const MatrixFree<dim, number> &
  get_system_matrix_free() const;

  const AlignedVector<VectorizedArray<number>>
  get_element_size() const;

  void
  compute_inverse_diagonal(VectorType &diagonal) const;

  void
  evaluate_non_linear_term(const VectorType &newton_step);

protected:
  virtual void
  do_cell_integral_local(FECellIntegrator &integrator) const;

  void
  do_cell_integral_range(
    const MatrixFree<dim, number> &              matrix_free,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &range) const;

  static IndexSet
  get_refinement_edges(const MatrixFree<dim, number> &matrix_free);

  MatrixFree<dim, number>                matrix_free;
  AffineConstraints<number>              constraints;
  mutable TrilinosWrappers::SparseMatrix system_matrix;
  SimulationParameters<dim>              parameters;
  AlignedVector<VectorizedArray<number>> element_size;
  Table<2, Tensor<1, dim + 1, VectorizedArray<number>>>
    nonlinear_previous_values;
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_gradient;
  Table<2, Tensor<1, dim + 1, Tensor<1, dim, VectorizedArray<number>>>>
    nonlinear_previous_hessian_diagonal;

  // Variables needed for local smoothing
  std::vector<unsigned int>                      constrained_indices;
  mutable std::vector<std::pair<number, number>> constrained_values;
  std::vector<unsigned int>                      edge_constrained_indices;
  bool has_edge_constrained_indices = false;
  mutable std::vector<std::pair<number, number>> edge_constrained_values;
  std::vector<bool>                              edge_constrained_cell;
};


template <int dim, typename number>
class NavierStokesSUPGPSPGOperator
  : public NavierStokesOperatorBase<dim, number>
{
public:
  using FECellIntegrator = FEEvaluation<dim, -1, 0, dim + 1, number>;
  NavierStokesSUPGPSPGOperator();

private:
  void
  do_cell_integral_local(FECellIntegrator &integrator) const override;
};

#endif
