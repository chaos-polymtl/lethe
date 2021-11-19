#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/vof_scratch_data.h>

template <int dim>
void
VOFScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_fs.get_quadrature().size();
  this->n_dofs     = fe_values_fs.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);


  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_gradient_values = std::vector<Tensor<2, dim>>(n_q_points);


  this->present_phase_values = std::vector<double>(n_q_points);
  this->phase_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->phase_laplacians     = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));



  // Velocity for SDIRK schemes
  this->stages_phase_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  this->phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
}


template class VOFScratchData<2>;
template class VOFScratchData<3>;
