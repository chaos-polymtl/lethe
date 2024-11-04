// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_phase_gradient_projection.h>

template <int dim>
void
VOFPhaseGradientProjection<dim>::assemble_system_matrix_and_rhs()
{
  // Reinitialize system matrix and right-hand side (rhs)
  this->system_matrix = 0;
  this->system_rhs    = 0;

  // Get VOF DoFHandler
  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics_interface->get_dof_handler(PhysicsID::VOF);

  // Initialize FEValues for phase fraction gradient projection and VOF
  FEValues<dim> fe_values_phase_gradient_projection(*this->mapping,
                                                    *this->fe,
                                                    *this->cell_quadrature,
                                                    update_values |
                                                      update_gradients |
                                                      update_JxW_values);
  FEValues<dim> fe_values_vof(*this->mapping,
                              dof_handler_vof->get_fe(),
                              *this->cell_quadrature,
                              update_gradients);

  // Initialize size of arrays
  const double n_q_points =
    fe_values_phase_gradient_projection.get_quadrature().size();
  const double n_dofs_per_cell =
    fe_values_phase_gradient_projection.get_fe().n_dofs_per_cell();

  // Initialize local matrix and rhs
  FullMatrix<double> local_matrix(n_dofs_per_cell, n_dofs_per_cell);
  Vector<double>     local_rhs(n_dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices(n_dofs_per_cell);

  // Extractor for phase fraction gradient vector
  FEValuesExtractors::Vector phase_fraction_gradients(0);

  // Initialize solution arrays
  std::vector<Tensor<1, dim>> present_phase_gradient_projection_values(
    n_q_points);
  std::vector<Tensor<2, dim>> present_phase_gradient_projection_gradients(
    n_q_points);
  std::vector<Tensor<1, dim>> present_filtered_vof_phase_gradients(n_q_points);

  // Initialize shape function arrays
  std::vector<Tensor<1, dim>> phi(n_dofs_per_cell);
  std::vector<Tensor<2, dim>> grad_phi(n_dofs_per_cell);

  // Get the diffusion factor
  const double diffusion_factor =
    this->simulation_parameters.multiphysics.vof_parameters
      .surface_tension_force.phase_fraction_gradient_diffusion_factor;

  // Loop over phase gradient projection cells
  for (const auto &cell : this->dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize local linear system
          local_matrix = 0;
          local_rhs    = 0;

          // Get VOF cell iterator
          typename DoFHandler<dim>::active_cell_iterator vof_cell(
            &(*this->triangulation),
            cell->level(),
            cell->index(),
            dof_handler_vof);

          // Reinitialize FEValues with corresponding cells
          fe_values_phase_gradient_projection.reinit(cell);
          fe_values_vof.reinit(vof_cell);

          // Get vector of Jacobi determinant times the quadrature weights
          std::vector<double> JxW_vec =
            fe_values_phase_gradient_projection.get_JxW_values();

          // Compute cell size
          auto &fe_phase_gradient_projection =
            fe_values_phase_gradient_projection.get_fe();
          const double h =
            compute_cell_diameter<dim>(compute_cell_measure_with_JxW(JxW_vec),
                                       fe_phase_gradient_projection.degree);

          // Get projected phase fraction gradient values, gradients and VOF
          // phase fraction gradients
          fe_values_phase_gradient_projection[phase_fraction_gradients]
            .get_function_values(this->present_solution,
                                 present_phase_gradient_projection_values);
          fe_values_phase_gradient_projection[phase_fraction_gradients]
            .get_function_gradients(
              this->present_solution,
              present_phase_gradient_projection_gradients);
          fe_values_vof.get_function_gradients(
            *this->multiphysics_interface->get_filtered_solution(
              PhysicsID::VOF),
            present_filtered_vof_phase_gradients);

          // Loop over quadrature points
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              // Shape functions
              for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
                {
                  phi[k] = fe_values_phase_gradient_projection
                             [phase_fraction_gradients]
                               .value(k, q);
                  grad_phi[k] = fe_values_phase_gradient_projection
                                  [phase_fraction_gradients]
                                    .gradient(k, q);
                }

              // Local linear system assembly
              for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
                {
                  // Assemble local matrix
                  for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                    {
                      local_matrix(i, j) +=
                        (phi[i] * phi[j] +
                         h * h * diffusion_factor *
                           scalar_product(grad_phi[i], grad_phi[j])) *
                        JxW_vec[q];
                    }
                  // Assemble local right-hand side (rhs)
                  local_rhs(i) += phi[i] *
                                  present_filtered_vof_phase_gradients[q] *
                                  JxW_vec[q];
                }
            }

          // Distribute the local contributions to the global system
          cell->get_dof_indices(local_dof_indices);
          this->constraints.distribute_local_to_global(local_matrix,
                                                       local_rhs,
                                                       local_dof_indices,
                                                       this->system_matrix,
                                                       this->system_rhs);
        }
    }
  this->system_matrix.compress(VectorOperation::add);
  this->system_rhs.compress(VectorOperation::add);
}


template class VOFPhaseGradientProjection<2>;
template class VOFPhaseGradientProjection<3>;
