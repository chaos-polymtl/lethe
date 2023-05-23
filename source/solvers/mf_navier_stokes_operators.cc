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

#include "solvers/mf_navier_stokes_operators.h"

template <int dim>
NavierStokesOperator<dim>::NavierStokesOperator()
  : MatrixFreeOperators::Base<dim, VectorType>()
{}

template <int dim>
void
NavierStokesOperator<dim>::clear()
{
  nonlinear_values.reinit(0, 0);
  MatrixFreeOperators::Base<dim, VectorType>::clear();
}

template <int dim>
void
NavierStokesOperator<dim>::reinit_operator_parameters()
{}

template <int dim>
void
NavierStokesOperator<dim>::evaluate_newton_step(const VectorType &newton_step)
{
  const unsigned int n_cells = this->data->n_cell_batches();
  FECellIntegrator   phi(*this->data);

  nonlinear_values.reinit(n_cells, phi.n_q_points);

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values_plain(newton_step);
      phi.evaluate(EvaluationFlags::values);

      for (unsigned int q = 0; q < phi.n_q_points; ++q)
        {
          // TODO: Evaluate non-linearity
          nonlinear_values(cell, q) = VectorizedArray<double>(0.0);
        }
    }
}

template <int dim>
void
NavierStokesOperator<dim>::compute_diagonal()
{
  this->inverse_diagonal_entries.reset(
    new DiagonalMatrix<LinearAlgebra::distributed::Vector<double>>());
  LinearAlgebra::distributed::Vector<double> &inverse_diagonal =
    this->inverse_diagonal_entries->get_vector();
  this->data->initialize_dof_vector(inverse_diagonal);

  MatrixFreeTools::compute_diagonal(*this->data,
                                    inverse_diagonal,
                                    &NavierStokesOperator::local_compute,
                                    this);

  for (auto &diagonal_element : inverse_diagonal)
    {
      diagonal_element =
        (std::abs(diagonal_element) > 1.0e-10) ? (1.0 / diagonal_element) : 1.0;
    }
}

template <int dim>
void
NavierStokesOperator<dim>::apply_add(VectorType &      dst,
                                     const VectorType &src) const
{
  this->data->cell_loop(&NavierStokesOperator::local_apply, this, dst, src);
}

template <int dim>
void
NavierStokesOperator<dim>::Tapply_add(VectorType &      dst,
                                      const VectorType &src) const
{
  this->data->cell_loop(&NavierStokesOperator::local_apply_transpose,
                        this,
                        dst,
                        src);
}

template <int dim>
void
NavierStokesOperator<dim>::local_apply(
  const MatrixFree<dim, double> &data,
  VectorType &                   dst,
  const VectorType & /* src */,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FECellIntegrator phi(data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0),
                      phi.get_matrix_free().n_cell_batches());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);

      phi.reinit(cell);

      // TODO: implement operator

      phi.integrate_scatter(EvaluationFlags::values |
                              EvaluationFlags::gradients,
                            dst);
    }
}

template <int dim>
void
NavierStokesOperator<dim>::local_apply_transpose(
  const MatrixFree<dim, double> &data,
  VectorType &                   dst,
  const VectorType & /* src */,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  FECellIntegrator phi(data);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      AssertDimension(nonlinear_values.size(0),
                      phi.get_matrix_free().n_cell_batches());
      AssertDimension(nonlinear_values.size(1), phi.n_q_points);

      phi.reinit(cell);

      // TODO: implement operator

      phi.integrate_scatter(EvaluationFlags::values |
                              EvaluationFlags::gradients,
                            dst);
    }
}

template <int dim>
void
NavierStokesOperator<dim>::local_compute(FECellIntegrator &phi) const
{
  AssertDimension(nonlinear_values.size(0),
                  phi.get_matrix_free().n_cell_batches());
  AssertDimension(nonlinear_values.size(1), phi.n_q_points);

  for (unsigned int q = 0; q < phi.n_q_points; ++q)
    {
      // TODO: implement
    }

  phi.integrate(EvaluationFlags::values | EvaluationFlags::gradients);
}

template class NavierStokesOperator<2>;
template class NavierStokesOperator<3>;