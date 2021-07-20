#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/navier_stokes_scratch_data.h>

template <int dim>
void
NavierStokesScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values.get_quadrature().size();
  this->n_dofs     = fe_values.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->rhs_force =
    std::vector<Vector<double>>(n_q_points, Vector<double>(dim + 1));
  this->force = std::vector<Tensor<1, dim>>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  this->pressure.component                = dim;
  // Velocity
  this->velocity_values      = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);
  this->velocity_gradients   = std::vector<Tensor<2, dim>>(n_q_points);
  this->velocity_laplacians  = std::vector<Tensor<1, dim>>(n_q_points);

  // Velocity for BDF schemes
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(n_q_points));

  // Velocity for SDIRK schemes
  this->stages_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    max_number_of_intermediary_stages(),
    std::vector<Tensor<1, dim>>(n_q_points));


  // Pressure
  this->pressure_values    = std::vector<double>(n_q_points);
  this->pressure_gradients = std::vector<Tensor<1, dim>>(n_q_points);

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_u = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->grad_phi_u = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->div_phi_u =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->hess_phi_u = std::vector<std::vector<Tensor<3, dim>>>(
    n_q_points, std::vector<Tensor<3, dim>>(n_dofs));
  this->laplacian_phi_u = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));

  // Pressure shape functions
  this->phi_p =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_p = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
}


template class NavierStokesScratchData<2>;
template class NavierStokesScratchData<3>;
