// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/fine_scale_velocity_history.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <string>

template <int dim>
void
FineScaleVelocityHistory<dim>::reinit(const Triangulation<dim> &triangulation,
                                      const unsigned int n_history_entries,
                                      const unsigned int n_q_points)
{
  this->n_history_entries = n_history_entries;
  this->n_q_points        = n_q_points;
  this->recording         = false;

  this->cell_data.assign(triangulation.n_active_cells(), {});
  this->new_values.assign(triangulation.n_active_cells(), {});
  if (this->n_values_per_cell() == 0)
    return;

  for (const auto &cell : triangulation.active_cell_iterators())
    if (cell->is_locally_owned())
      {
        this->cell_data[cell->active_cell_index()].assign(
          this->n_values_per_cell(), 0.);
        this->new_values[cell->active_cell_index()].assign(n_q_points * dim,
                                                           0.);
      }
}

template <int dim>
void
FineScaleVelocityHistory<dim>::advance()
{
  const unsigned int n_values_per_entry = this->n_q_points * dim;
  for (unsigned int cell = 0; cell < this->cell_data.size(); ++cell)
    {
      std::vector<double> &values = this->cell_data[cell];
      if (values.empty())
        continue;

      // Shift the entries by one (e.g. [2] <- [1] <- [0] for BDF3) and store
      // the new fine-scale velocity in entry 0.
      std::copy_backward(values.begin(),
                         values.end() - n_values_per_entry,
                         values.end());
      std::copy(this->new_values[cell].begin(),
                this->new_values[cell].end(),
                values.begin());
    }
}

template <int dim>
void
FineScaleVelocityHistory<dim>::set_cell_data(const CellData &new_cell_data)
{
  AssertThrow(new_cell_data.size() == this->cell_data.size(),
              ExcMessage(
                "The cell-wise data of the history of the dynamic RBVMS "
                "fine-scale velocity holds " +
                std::to_string(new_cell_data.size()) +
                " cells, but the history was allocated for " +
                std::to_string(this->cell_data.size()) + " active cells."));

  for (unsigned int cell = 0; cell < this->cell_data.size(); ++cell)
    if (!this->cell_data[cell].empty())
      {
        AssertDimension(new_cell_data[cell].size(),
                        this->cell_data[cell].size());
        this->cell_data[cell] = new_cell_data[cell];
      }
}

template class FineScaleVelocityHistory<2>;
template class FineScaleVelocityHistory<3>;
