#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/heat_transfer_scratch_data.h>

template <int dim>
void
HeatTransferScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_T.get_quadrature().size();
  this->n_dofs     = fe_values_T.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_gradient_values = std::vector<Tensor<2, dim>>(n_q_points);
  // Temperature
  this->phi_T           = std::vector<double>(n_dofs);
  this->grad_phi_T      = std::vector<Tensor<1, dim>>(n_dofs);
  this->hess_phi_T      = std::vector<Tensor<2, dim>>(n_dofs);
  this->laplacian_phi_T = std::vector<double>(n_dofs);

  this->present_temperature_values = std::vector<double>(n_q_points);
  this->temperature_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->present_temperature_laplacians = std::vector<double>(n_q_points);

  this->previous_T_values = std::vector<std::vector<double>>(n_q_points);
  this->previous_T_gradients =
    std::vector<std::vector<Tensor<1, dim>>>(n_q_points);


  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  this->phase_values = std::vector<double>(n_q_points);
}


template class HeatTransferScratchData<2>;
template class HeatTransferScratchData<3>;
