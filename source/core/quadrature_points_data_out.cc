// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/quadrature_points_data_out.h>

#include <deal.II/fe/fe_values.h>

template <int dim>
void
QuadraturePointsDataOut<dim>::build_patches(
  const dealii::DoFHandler<dim> &dof_handler,
  const dealii::Mapping<dim>    &mapping,
  const dealii::Quadrature<dim> &quadrature,
  const std::function<
    double(const typename dealii::DoFHandler<dim>::active_cell_iterator &)>
    &filter_radius)
{
  patches.clear();
  // Upper bound: also counts the ghost and artificial cells stored locally
  patches.reserve(dof_handler.get_triangulation().n_active_cells() *
                  quadrature.size());

  dealii::FEValues<dim> fe_values(mapping,
                                  dof_handler.get_fe(),
                                  quadrature,
                                  dealii::update_quadrature_points);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);
      const double radius = filter_radius(cell);

      for (const auto &q : fe_values.get_quadrature_points())
        {
          dealii::DataOutBase::Patch<0, dim> patch;
          patch.vertices[0] = q;
          patch.data.reinit(dataset_names.size(), 1);
          patch.data(0, 0) = radius;
          patch.data(1, 0) = 2. * radius;
          patch.data(2, 0) = cell->global_active_cell_index();
          patches.push_back(patch);
        }
    }
}

template class QuadraturePointsDataOut<2>;
template class QuadraturePointsDataOut<3>;
