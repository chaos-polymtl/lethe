#include <solvers/vof_algebraic_reinitialization_scratch_data.h>

template <int dim>
void
VOFAlgebraicReinitializationScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_vof.get_quadrature().size();
  this->n_dofs     = fe_values_vof.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(this->n_q_points);

  // Initialize solution arrays
  this->present_phase_algebraic_reinit_values =
    std::vector<double>(this->n_q_points);
  this->present_phase_algebraic_reinit_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_phase_algebraic_reinit_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));

  // Initialize arrays related to VOF phase gradient for normal vector
  // computation
  this->present_vof_phase_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Initialize arrays related to shape functions
  this->ksi =
    std::vector<std::vector<double>>(this->n_q_points,
                                     std::vector<double>(this->n_dofs));
  this->grad_ksi = std::vector<std::vector<Tensor<1, dim>>>(
    this->n_q_points, std::vector<Tensor<1, dim>>(this->n_dofs));
}

template class VOFAlgebraicReinitializationScratchData<2>;
template class VOFAlgebraicReinitializationScratchData<3>;
