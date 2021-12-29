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

  // Initialize arrays related to temperature
  this->present_temperature_values = std::vector<double>(n_q_points);
  this->temperature_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->present_temperature_laplacians = std::vector<double>(n_q_points);

  // Initial arrays for physical properties
  this->specific_heat        = std::vector<double>(n_q_points);
  this->thermal_conductivity = std::vector<double>(n_q_points);
  this->density              = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_temperature_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  this->previous_temperature_gradients =
    std::vector<std::vector<Tensor<1, dim>>>(
      maximum_number_of_previous_solutions(),
      std::vector<Tensor<1, dim>>(n_q_points));

  // Velocity for SDIRK schemes
  this->stages_temperature_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_T =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_T = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_T = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_T =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  // Allocate memory for the physical properties
  specific_heat        = std::vector<double>(n_q_points);
  density              = std::vector<double>(n_q_points);
  thermal_conductivity = std::vector<double>(n_q_points);
}

template <int dim>
void
HeatTransferScratchData<dim>::enable_vof(const FiniteElement<dim> &fe,
                                         const Quadrature<dim> &   quadrature,
                                         const Mapping<dim> &      mapping)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // VOF
  phase_values = std::vector<double>(this->n_q_points);
  previous_vof_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
}


template class HeatTransferScratchData<2>;
template class HeatTransferScratchData<3>;
