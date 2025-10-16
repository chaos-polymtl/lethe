// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>

#include <solvers/vof_scratch_data.h>

template <int dim>
void
VOFScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_vof.get_quadrature().size();
  this->n_dofs     = fe_values_vof.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities_fd.first_vector_component = 0;
  // Velocity
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  this->velocity_gradient_values =
    std::vector<Tensor<2, dim>>(this->n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);

  this->present_phase_values = std::vector<double>(this->n_q_points);
  this->phase_gradients      = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_phase_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  this->phase_laplacians = std::vector<double>(this->n_q_points);

  // Velocity for BDF schemes
  this->previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));

  // Initialize arrays related to shape functions
  this->phi.reinit(n_q_points, n_dofs);
  this->grad_phi.reinit(n_q_points, n_dofs);
  this->hess_phi.reinit(n_q_points, n_dofs);
  this->laplacian_phi.reinit(n_q_points, n_dofs);
}

template class VOFScratchData<2>;
template class VOFScratchData<3>;
