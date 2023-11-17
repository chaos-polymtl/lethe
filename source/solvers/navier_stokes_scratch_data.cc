#include <core/bdf.h>
#include <core/dem_properties.h>
#include <core/utilities.h>

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

  // Initialize component array
  this->components = std::vector<unsigned int>(n_dofs);

  // Forcing term array
  this->rhs_force =
    std::vector<Vector<double>>(n_q_points, Vector<double>(dim + 1));
  this->force       = std::vector<Tensor<1, dim>>(n_q_points);
  this->mass_source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  this->pressure.component                = dim;

  // Velocity
  this->velocity_values      = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);
  this->velocity_gradients   = std::vector<Tensor<2, dim>>(n_q_points);
  this->velocity_laplacians  = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_hessians    = std::vector<Tensor<3, dim>>(n_q_points);
  this->shear_rate           = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(n_q_points));

  // Pressure
  this->pressure_values         = std::vector<double>(n_q_points);
  this->pressure_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->pressure_scaling_factor = 1;

  // Pressure for BDF schemes (compressible Navier-Stokes)
  this->previous_pressure_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
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

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));

  density                             = std::vector<double>(n_q_points);
  dynamic_viscosity                   = std::vector<double>(n_q_points);
  kinematic_viscosity                 = std::vector<double>(n_q_points);
  thermal_expansion                   = std::vector<double>(n_q_points);
  grad_kinematic_viscosity_shear_rate = std::vector<double>(n_q_points);

  previous_density =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_vof(
  const FiniteElement<dim>          &fe,
  const Quadrature<dim>             &quadrature,
  const Mapping<dim>                &mapping,
  const Parameters::VOF_PhaseFilter &phase_filter_parameters)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate VOF values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                  = std::vector<double>(n_q_points);
  density_1                  = std::vector<double>(n_q_points);
  dynamic_viscosity_0        = std::vector<double>(n_q_points);
  dynamic_viscosity_1        = std::vector<double>(n_q_points);
  thermal_expansion_0        = std::vector<double>(n_q_points);
  thermal_expansion_1        = std::vector<double>(n_q_points);
  surface_tension            = std::vector<double>(n_q_points);
  surface_tension_gradient   = std::vector<double>(n_q_points);
  compressibility_multiplier = std::vector<double>(n_q_points);

  // Create filter
  filter = VolumeOfFluidFilterBase::model_cast(phase_filter_parameters);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_vof(
  const FiniteElement<dim>                       &fe,
  const Quadrature<dim>                          &quadrature,
  const Mapping<dim>                             &mapping,
  const std::shared_ptr<VolumeOfFluidFilterBase> &filter)
{
  gather_vof    = true;
  fe_values_vof = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate VOF values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                  = std::vector<double>(n_q_points);
  density_1                  = std::vector<double>(n_q_points);
  dynamic_viscosity_0        = std::vector<double>(n_q_points);
  dynamic_viscosity_1        = std::vector<double>(n_q_points);
  thermal_expansion_0        = std::vector<double>(n_q_points);
  thermal_expansion_1        = std::vector<double>(n_q_points);
  surface_tension            = std::vector<double>(n_q_points);
  surface_tension_gradient   = std::vector<double>(n_q_points);
  compressibility_multiplier = std::vector<double>(n_q_points);

  // Create filter
  this->filter = filter;
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_cahn_hilliard(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_cahn_hilliard    = true;
  fe_values_cahn_hilliard = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CahnHilliard values
  phase_order_cahn_hilliard_values = std::vector<double>(this->n_q_points);
  chemical_potential_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  // Allocate CahnHilliard gradients
  phase_order_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  chemical_potential_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));

  // Allocate physical properties
  density_0                = std::vector<double>(n_q_points);
  density_1                = std::vector<double>(n_q_points);
  dynamic_viscosity_0      = std::vector<double>(n_q_points);
  dynamic_viscosity_1      = std::vector<double>(n_q_points);
  thermal_expansion_0      = std::vector<double>(n_q_points);
  thermal_expansion_1      = std::vector<double>(n_q_points);
  surface_tension          = std::vector<double>(n_q_points);
  surface_tension_gradient = std::vector<double>(n_q_points);
  // mobility_cahn_hilliard   = std::vector<double>(n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_projected_phase_fraction_gradient(
  const FiniteElement<dim> &fe_projected_phase_fraction_gradient,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_projected_phase_fraction_gradient = true;
  fe_values_projected_phase_fraction_gradient =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_projected_phase_fraction_gradient,
                                    quadrature,
                                    update_values | update_gradients);

  // phase fraction gradient (PFG)
  projected_phase_fraction_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_curvature(
  const FiniteElement<dim> &fe_curvature,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_curvature    = true;
  fe_values_curvature = std::make_shared<FEValues<dim>>(
    mapping, fe_curvature, quadrature, update_values | update_gradients);

  // curvature
  curvature_values = std::vector<double>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_void_fraction(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_void_fraction    = true;
  fe_values_void_fraction = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Void Fraction
  void_fraction_values = std::vector<double>(this->n_q_points);
  previous_void_fraction_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  void_fraction_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_particle_fluid_interactions(
  const unsigned int n_global_max_particles_per_cell,
  const bool         enable_void_fraction_interpolation)
{
  gather_particles_information     = true;
  max_number_of_particles_per_cell = n_global_max_particles_per_cell;
  interpolated_void_fraction       = enable_void_fraction_interpolation;

  // Velocities
  particle_velocity =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  fluid_velocity_at_particle_location =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  cell_void_fraction = std::vector<double>(n_global_max_particles_per_cell);
  fluid_particle_relative_velocity_at_particle_location =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  Re_particle = std::vector<double>(n_global_max_particles_per_cell);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_heat_transfer(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_temperature    = true;
  fe_values_temperature = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  temperature_values    = std::vector<double>(this->n_q_points);
  temperature_gradients = std::vector<Tensor<1, dim>>(this->n_q_points);
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
}


template <int dim>
void
NavierStokesScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::temperature) &&
      gather_temperature)
    {
      set_field_vector(field::temperature,
                       this->temperature_values,
                       this->fields);
    }

  if (properties_manager.field_is_required(field::pressure))
    {
      set_field_vector(field::pressure, this->pressure_values, this->fields);
    }

  if (properties_manager.field_is_required(field::shear_rate))
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Calculate shear rate (at each q)
          const Tensor<2, dim> shear_rate_tensor =
            velocity_gradients[q] + transpose(velocity_gradients[q]);

          // Calculate the shear rate magnitude
          shear_rate[q] = calculate_shear_rate_magnitude(shear_rate_tensor);
        }

      set_field_vector(field::shear_rate, shear_rate, this->fields);
    }

  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          // In the case of an incompressible flow, viscosity is the only
          // required property
          const auto rheology_model = properties_manager.get_rheology();
          rheology_model->vector_value(fields, kinematic_viscosity);
          kinematic_viscosity_scale =
            rheology_model->get_kinematic_viscosity_scale();

          if (!properties_manager.density_is_constant())
            {
              // For a weakly compressible flow, density variations will play a
              // role
              const auto density_model = properties_manager.get_density();
              density_model->vector_value(fields, density);
              density_psi = density_model->get_psi();
              density_ref = density_model->get_density_ref();
              // Dynamic viscosity is also necessary for compressible flows
              rheology_model->get_dynamic_viscosity_vector(density_ref,
                                                           fields,
                                                           dynamic_viscosity);
            }

          if (properties_manager.is_non_newtonian())
            {
              // Calculate derivative of viscosity with respect to shear rate
              rheology_model->vector_jacobian(
                fields, field::shear_rate, grad_kinematic_viscosity_shear_rate);
            }

          if (gather_temperature)
            {
              const auto thermal_expansion_model =
                properties_manager.get_thermal_expansion();
              thermal_expansion_model->vector_value(fields, thermal_expansion);
            }
          break;
        }
      case 2:
        {
          // We need both density and viscosity
          const auto density_model_0  = properties_manager.get_density(0);
          const auto rheology_model_0 = properties_manager.get_rheology(0);
          const auto density_model_1  = properties_manager.get_density(1);
          const auto rheology_model_1 = properties_manager.get_rheology(1);

          density_ref_0 = density_model_0->get_density_ref();
          density_ref_1 = density_model_1->get_density_ref();

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
              // Gather surface tension gradient only if necessary
              if (!properties_manager.surface_tension_is_constant())
                surface_tension_model->vector_jacobian(
                  fields, field::temperature, surface_tension_gradient);
            }

          density_model_0->vector_value(fields, density_0);
          rheology_model_0->get_dynamic_viscosity_vector(density_ref_0,
                                                         fields,
                                                         dynamic_viscosity_0);

          density_model_1->vector_value(fields, density_1);
          rheology_model_1->get_dynamic_viscosity_vector(density_ref_1,
                                                         fields,
                                                         dynamic_viscosity_1);

          if (gather_temperature)
            {
              const auto thermal_expansion_model_0 =
                properties_manager.get_thermal_expansion(0);
              const auto thermal_expansion_model_1 =
                properties_manager.get_thermal_expansion(1);
              thermal_expansion_model_0->vector_value(fields,
                                                      thermal_expansion_0);
              thermal_expansion_model_1->vector_value(fields,
                                                      thermal_expansion_1);
            }

          if (gather_vof && !gather_cahn_hilliard)
            {
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  double filtered_phase_value = this->filtered_phase_values[q];

                  density[q] = calculate_point_property(filtered_phase_value,
                                                        this->density_0[q],
                                                        this->density_1[q]);

                  dynamic_viscosity[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->dynamic_viscosity_0[q],
                                             this->dynamic_viscosity_1[q]);

                  thermal_expansion[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->thermal_expansion_0[q],
                                             this->thermal_expansion_1[q]);
                }

              // Gather density_psi for isothermal compressible NS equations
              if (!properties_manager.density_is_constant())
                {
                  density_psi_0 = density_model_0->get_psi();
                  density_psi_1 = density_model_1->get_psi();

                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      double filtered_phase_value =
                        this->filtered_phase_values[q];
                      // To avoid calculating twice in the assemblers
                      compressibility_multiplier[q] =
                        filtered_phase_value *
                          (density_0[q] / density_1[q] * density_psi_1 -
                           density_psi_0) +
                        density_psi_0;
                    }
                }

              for (unsigned p = 0; p < previous_phase_values.size(); ++p)
                {
                  // Blend the physical properties using the VOF field
                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      // Calculate previous density (right now assumes constant
                      // density model per phase)
                      previous_density[p][q] = calculate_point_property(
                        filter->filter_phase(this->previous_phase_values[p][q]),
                        this->density_0[q],
                        this->density_1[q]);
                    }
                }
              break;
            }
          else if (gather_cahn_hilliard && !gather_vof)
            {
              // std::cout<< "Computing properties for Cahn-Hilliard" <<
              // std::endl;
              //  Blend the physical properties using the CahnHilliard field
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  // this->density_diff =
                  //   0.5 * std::abs(density_0[q] - density_1[q]);

                  double phase_order_cahn_hilliard_value =
                    this->phase_order_cahn_hilliard_values[q];

                  /*
                  std::cout<<"phase value = "
                  <<phase_order_cahn_hilliard_value<<std::endl;
                  std::cout<<"density 0 = " <<density_0[q]<<std::endl;
                  std::cout<<"density 1 = " <<density_1[q]<<std::endl;
                  std::cout<<"expected density at point q = "<<
                  0.5*density_0[q]*(phase_order_cahn_hilliard_value+1)
                  -0.5*density_1[q]*(phase_order_cahn_hilliard_value-1)<<std::endl;
*/
                  density[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->density_0[q],
                    this->density_1[q]);

                  /*
                  std::cout<<"computed density at point q = "<<
                  density[q]<<std::endl;

                  std::cout<<"mu 0 = " <<dynamic_viscosity_0[q]<<std::endl;
                  std::cout<<"mu  1 = " <<dynamic_viscosity_1[q]<<std::endl;
                  std::cout<<"expected mu at point q = "<<
                  0.5*dynamic_viscosity_0[q]*(phase_order_cahn_hilliard_value+1)
                  -0.5*dynamic_viscosity_1[q]*(phase_order_cahn_hilliard_value-1)<<std::endl;
*/
                  dynamic_viscosity[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->dynamic_viscosity_0[q],
                    this->dynamic_viscosity_1[q]);

                  // std::cout<<"computed mu at point q = "<<
                  // dynamic_viscosity[q]<<std::endl;
                  /*
                const auto material_interaction_id =
                  properties_manager.get_material_interaction_id(
                    material_interactions_type::fluid_fluid, 0, 1);
                const auto mobility_cahn_hilliard_model =
                  properties_manager.get_mobility_cahn_hilliard(
                    material_interaction_id);
                mobility_cahn_hilliard_model->vector_value(
                  fields, mobility_cahn_hilliard);
                  std::cout<<"mobility = "
                <<mobility_cahn_hilliard[q]<<std::endl;
                   */
                }
              break;
            }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template class NavierStokesScratchData<2>;
template class NavierStokesScratchData<3>;
