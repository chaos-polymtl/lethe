// SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/lethe_grid_tools.h>

#include <solvers/stabilization.h>

template <int dim>
void
moe_scalar_limiter(const DoFHandler<dim> &dof_handler,
                   GlobalVectorType      &locally_relevant_vector,
                   GlobalVectorType      &locally_owned_vector)
{
  // The limiter requires a vector in which the value can be modified.
  // For this we use the local_evaluation_point. We thus begin by making sure
  // that the local_evaluation_point is set at the present value of the solution
  // to limit.
  locally_owned_vector = locally_relevant_vector;

  // We define the cutoff function from the Moe paper using a lambda function.
  // It's a bit weird this cutoff function, but YOLO I guess. Why do we divide y
  // by 1.1? I do not know, I believe it to be one of the mysteries of life.
  auto phi_limit = [](const double y) { return std::min(1., y / 1.1); };

  // The strategy for limiting requires looping over the cells twice.
  // First loop over all active cells (local and ghost):
  // 1. Calculate the max and the min value of the field we wish to limit
  // Loop over every local cells:
  // 2. Calculate approximate upper and lower bounds using the neighbors
  // 3. Calculate the value of the theta limiting parameter
  // 4. Rescale the nodal values of the function using the average value of the
  // solution within the element and the theta limiter.

  // We need a vertices to cell map to have access to the neighbors rapidly
  // We obtain this map using the LetheGridTools functionnalities.
  std::map<unsigned int,
           std::set<typename DoFHandler<dim>::active_cell_iterator>>
    vertices_to_cell;
  LetheGridTools::vertices_cell_mapping(dof_handler, vertices_to_cell);

  // Step 1, loop over every cell and calculate the max, the min and the mean.
  // These values are stored in std::map for rapid access.
  std::map<typename ::dealii::types::global_cell_index, double>
    max_value_per_cells;
  std::map<typename ::dealii::types::global_cell_index, double>
    min_value_per_cells;
  std::map<typename ::dealii::types::global_cell_index, double>
    mean_value_per_cells;

  std::vector<types::global_dof_index> local_dof_indices(
    dof_handler.get_fe().n_dofs_per_cell());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // We need to loop over cells that are either locally_owned or ghost to
      // gather all the neighbors.
      if (cell->is_locally_owned() || cell->is_ghost())
        {
          cell->get_dof_indices(local_dof_indices);
          double max_value  = -DBL_MAX;
          double min_value  = DBL_MAX;
          double mean_value = 0;
          for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
            {
              // Get the max and the min value of the solution
              double dof_value = locally_relevant_vector(local_dof_indices[i]);
              max_value        = std::max(max_value, dof_value);
              min_value        = std::min(min_value, dof_value);
              // Calculate mean solution \bar{w}. Right now this is an average
              // of the nodal value which is lazy as balls. A better one will be
              // to do an integral, but this is dev and I am a lazy bum for now.
              mean_value += dof_value;
            }
          mean_value /= local_dof_indices.size();
          min_value_per_cells[cell->global_active_cell_index()]  = min_value;
          max_value_per_cells[cell->global_active_cell_index()]  = max_value;
          mean_value_per_cells[cell->global_active_cell_index()] = mean_value;
        }
    }

  // Loop over every cell and do:
  // 2. Calculate approximate upper and lower bound using the neighbhors
  // 3. Calculate the theta limiter
  // 4. Rescale the nodal values of the solution

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Get an approximate cell size
          // const double alpha = 50 * cell->diameter();
          cell->get_dof_indices(local_dof_indices);

          // Active neighbors include the current cell as well
          auto active_neighbors =
            LetheGridTools::find_cells_around_cell<dim>(vertices_to_cell, cell);


          // 2. Calculate approximate upper and lower bound using the neighbhors

          // Here we assume that alpha = 0 (see the Moe paper). The definition
          // of alpha in the original paper is flaky in my opinion, let's use
          // the more robust version as a started point.
          double upper_bound =
            mean_value_per_cells.at(cell->global_active_cell_index());
          double lower_bound =
            mean_value_per_cells.at(cell->global_active_cell_index());


          for (const auto &neighbor : active_neighbors)
            {
              // The shock capture mechanism theta calculation should not
              // consider the cell itself.
              if (cell->global_active_cell_index() ==
                  neighbor->global_active_cell_index())
                continue;

              // The original paper uses the max and the min of the neighbor
              // solutions. I have also tested with the mean values of the
              // neighbor. It seemed a bit more robust, but I'd rather follow
              // the paper for now. Let's just keep in mind that this can be a
              // fallback.

              upper_bound = std::max(upper_bound,
                                     max_value_per_cells.at(
                                       neighbor->global_active_cell_index()));

              lower_bound = std::min(lower_bound,
                                     min_value_per_cells.at(
                                       neighbor->global_active_cell_index()));
            }

          // 3. Calculate the value of the theta limiting using the max, min and
          // mean values as well as the bounds.
          const double max_value =
            max_value_per_cells.at(cell->global_active_cell_index());
          const double min_value =
            min_value_per_cells.at(cell->global_active_cell_index());
          const double mean_value =
            mean_value_per_cells.at(cell->global_active_cell_index());
          double y_max =
            (upper_bound - mean_value) / ((max_value - mean_value) + 1e-16);
          double y_min =
            (lower_bound - mean_value) / ((min_value - mean_value) + 1e-16);

          double theta = std::min(1., phi_limit(y_max));
          theta        = std::min(theta, phi_limit(y_min));

          // 4. Rescale the solution within the element
          for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
            {
              // Get the value of the solution at the DOFs and rescale it
              // using the mean value and the theta limiter.
              const double dof_value =
                locally_relevant_vector(local_dof_indices[i]);
              locally_owned_vector(local_dof_indices[i]) =
                mean_value + theta * (dof_value - mean_value);
            }
        }
    }

  // Reset the present solution and the evaluation point to the new updated
  // solution.
  locally_relevant_vector = locally_owned_vector;
}


template void
moe_scalar_limiter(const DoFHandler<2> &dof_handler,
                   GlobalVectorType    &locally_relevant_vector,
                   GlobalVectorType    &locally_owned_vector);

template void
moe_scalar_limiter(const DoFHandler<3> &dof_handler,
                   GlobalVectorType    &locally_relevant_vector,
                   GlobalVectorType    &locally_owned_vector);
