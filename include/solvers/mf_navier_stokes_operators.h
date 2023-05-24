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

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/tools.h>

using namespace dealii;

template <int dim>
class NavierStokesOperator
  : public MatrixFreeOperators::Base<dim,
                                     LinearAlgebra::distributed::Vector<double>>
{
public:
  using VectorType = LinearAlgebra::distributed::Vector<double>;

  using FECellIntegratorVelocity = FEEvaluation<dim, -1, 0, dim, double>;

  using FECellIntegratorPressure = FEEvaluation<dim, -1, 0, 1, double>;

  NavierStokesOperator();

  virtual void
  clear() override;

  void
  reinit_operator_parameters();

  void
  evaluate_newton_step(const VectorType &newton_step);

  virtual void
  compute_diagonal() override;

private:
  virtual void
  apply_add(VectorType &dst, const VectorType &src) const override;

  virtual void
  Tapply_add(VectorType &dst, const VectorType &src) const override;

  void
  local_apply(const MatrixFree<dim, double> &              data,
              VectorType &                                 dst,
              const VectorType &                           src,
              const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  local_apply_transpose(
    const MatrixFree<dim, double> &              data,
    VectorType &                                 dst,
    const VectorType &                           src,
    const std::pair<unsigned int, unsigned int> &cell_range) const;

  void
  local_compute(FECellIntegratorVelocity &integrator_v) const;

  Table<2, VectorizedArray<double>> nonlinear_values;
};


#endif