#include <core/bdf.h>

#include <solvers/cahn_hilliard_scratch_data.h>
template <int dim>
void
CahnHilliardScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_cahn_hilliard.get_quadrature().size();
  this->n_dofs     = fe_values_cahn_hilliard.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source_phase_order        = std::vector<double>(n_q_points);
  this->source_chemical_potential = std::vector<double>(n_q_points);

  // Initialize arrays related to phase order and chemical potential
  this->phase_order.component        = 0;
  this->chemical_potential.component = 1;

  // Variables for the Cahn-Hilliard equations
  this->phase_order_values            = std::vector<double>(n_q_points);
  this->phase_order_gradients         = std::vector<Tensor<1, dim>>(n_q_points);
  this->phase_order_laplacians        = std::vector<double>(n_q_points);
  this->chemical_potential_values     = std::vector<double>(n_q_points);
  this->chemical_potential_gradients  = std::vector<Tensor<1, dim>>(n_q_points);
  this->chemical_potential_laplacians = std::vector<double>(n_q_points);

  // Phase order for BDF schemes
  this->previous_phase_order_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  this->previous_chemical_potential_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Phase-order shape functions
  this->phi_phase =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_phase = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_phase = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_phase =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  // Chemical potential shape functions
  this->phi_potential =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));
  this->grad_phi_potential = std::vector<std::vector<Tensor<1, dim>>>(
    n_q_points, std::vector<Tensor<1, dim>>(n_dofs));
  this->hess_phi_potential = std::vector<std::vector<Tensor<2, dim>>>(
    n_q_points, std::vector<Tensor<2, dim>>(n_dofs));
  this->laplacian_phi_potential =
    std::vector<std::vector<double>>(n_q_points, std::vector<double>(n_dofs));

  // Velocity values
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  this->velocity_gradient_values =
    std::vector<Tensor<2, dim>>(this->n_q_points);

  // Surface tension values
    this->surface_tension = std::vector<double>(n_q_points);
}

template <int dim>
void
CahnHilliardScratchData<dim>::calculate_physical_properties()
{
    switch (properties_manager.get_number_of_fluids())
    {
        case 1:
        {
            throw std::runtime_error("Unsupported number of fluids (1)");
            break;
        }

        case 2:
        {
            // Gather properties from material interactions if necessary
            if (properties_manager.get_number_of_material_interactions() > 0)
            {
                const auto material_interaction_id =
                        properties_manager.get_material_interaction_id(
                                material_interactions_type::fluid_fluid, 0, 1);
                // Gather surface tension
                const auto surface_tension_model =
                        properties_manager.get_surface_tension(material_interaction_id);
                surface_tension_model->vector_value(fields, surface_tension);
            }
            break;
        }

        default:
            throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template class CahnHilliardScratchData<2>;
template class CahnHilliardScratchData<3>;
