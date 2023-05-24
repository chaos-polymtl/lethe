#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/cahn_hilliard_scratch_data.h>
template <int dim>
void
CahnHilliardScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_ch.get_quadrature().size();
  this->n_dofs     = fe_values_ch.get_fe().n_dofs_per_cell();

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
  // Velocity for SDIRK schemes
  this->stages_phase_order_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
                                     std::vector<double>(n_q_points));

  this->stages_chemical_potential_values =
    std::vector<std::vector<double>>(max_number_of_intermediary_stages(),
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
}

/*
 * DO NOT USE TWO FLUIDS WITH CAHN-HILLIARD FOR THE MOMENT, IT'LL ONLY PRODUCE A
 * SEGFAULT SINCE IT'LL TRY TO RERIEVE UN-INITIALIZED VARIABLES SINCE VOF IS
 * DESACTIVATED
 */
template <int dim>
void
CahnHilliardScratchData<dim>::calculate_physical_properties()
{
  if (material_id < 1 || properties_manager.get_number_of_solids() < 1)
    {
      switch (properties_manager.get_number_of_fluids())
        {
          case 1:
            {
              const auto density_model  = properties_manager.get_density();
              const auto rheology_model = properties_manager.get_rheology();

              density_model->vector_value(fields, density);
              rheology_model->vector_value(fields, viscosity);

              break;
            }
          case 2:
            {
              const auto density_models =
                properties_manager.get_density_vector();
              const auto rheology_models =
                properties_manager.get_rheology_vector();

              density_models[0]->vector_value(fields, density_0);
              rheology_models[0]->vector_value(fields, viscosity_0);

              density_models[1]->vector_value(fields, density_1);
              rheology_models[1]->vector_value(fields, viscosity_1);

              // Blend the physical properties using the phase field
              // (Cahn-Hilliard)
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  auto phase_order_value = this->phase_order_values[q];

                  density[q] = calculate_point_property_ch(phase_order_value,
                                                           this->density_0[q],
                                                           this->density_1[q]);

                  viscosity[q] =
                    calculate_point_property_ch(phase_order_value,
                                                this->viscosity_0[q],
                                                this->viscosity_1[q]);
                }
            }
            break;
          default:
            throw std::runtime_error("Unsupported number of fluids (>2)");
        }
    }
  else
    {
      const auto density_model = properties_manager.get_density(0, material_id);
      const auto rheology_model =
        properties_manager.get_rheology(0, material_id);

      density_model->vector_value(fields, density);
      rheology_model->vector_value(fields, viscosity);
    }
}

template class CahnHilliardScratchData<2>;
template class CahnHilliardScratchData<3>;
