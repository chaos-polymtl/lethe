// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/quadrature_points_data_out.h>

template <int dim>
void
QuadraturePointsDataOut<dim>::build_patches(
  const dealii::Point<dim>              &location,
  const double                           filter_radius,
  const dealii::types::global_cell_index cell_index)
{
  dealii::DataOutBase::Patch<0, dim> patch;
  patch.vertices[0] = location;
  patch.data.reinit(dataset_names.size(), 1);
  patch.data(0, 0) = filter_radius;
  patch.data(1, 0) = 2. * filter_radius;
  patch.data(2, 0) = cell_index;
  patches.push_back(patch);
}

template class QuadraturePointsDataOut<2>;
template class QuadraturePointsDataOut<3>;