#include <core/bdf.h>

#include <solvers/tracer_scratch_data.h>

template <int dim>
void
TracerScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_tracer.get_quadrature().size();
  this->n_dofs     = fe_values_tracer.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW = std::vector<double>(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values = std::vector<Tensor<1, dim>>(n_q_points);
  // Tracer
  this->tracer_values        = std::vector<double>(n_q_points);
  this->tracer_gradients     = std::vector<Tensor<1, dim>>(n_q_points);
  this->tracer_laplacians    = std::vector<double>(n_q_points);
  this->tracer_diffusivity   = std::vector<double>(n_q_points);
  this->tracer_diffusivity_0 = std::vector<double>(n_q_points);
  this->tracer_diffusivity_1 = std::vector<double>(n_q_points);

  // Solid signed distance function
  if (properties_manager.field_is_required(field::levelset))
    this->sdf_values = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_tracer_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

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

  if (properties_manager.field_is_required(field::levelset))
    fields.insert(
      std::pair<field, std::vector<double>>(field::levelset, n_q_points));
}


template <int dim>
void
TracerScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::levelset))
    set_field_vector(field::levelset, this->sdf_values, this->fields);

  switch (properties_manager.get_number_of_solids())
    {
      // No solid
      case 0:
        {
          switch (properties_manager.get_number_of_fluids())
            {
              // One fluid
              case 1:
                {
                  // In this case, only viscosity is the required property
                  const auto diffusivity_model =
                    properties_manager.get_tracer_diffusivity();
                  diffusivity_model->vector_value(fields, tracer_diffusivity);
                  break;
                }
                // Two fluids
              case 2:
                {
                  // In this case, we need both density and viscosity
                  const auto diffusivity_models =
                    properties_manager.get_tracer_diffusivity_vector();

                  diffusivity_models[0]->vector_value(fields,
                                                      tracer_diffusivity_0);
                  diffusivity_models[1]->vector_value(fields,
                                                      tracer_diffusivity_1);

                  // TODO Incomplete at the present time because the tracer VOF
                  // complete is not finished Blend the physical properties
                  // using the VOF field
                  // for (unsigned int q = 0; q < this->n_q_points; ++q)
                    // {
                      //          tracer_diffusivity[q] =
                      //            calculate_point_property(this->phase_values[q],
                      //                                     this->density_0[q],
                      //                                     this->density_1[q]);
                    // }
                  break;
                }
              default:
                throw std::runtime_error("Unsupported number of fluids (>2)");
            }
          break;
        }
        // One solid
      case 1:
        {
          switch (properties_manager.get_number_of_fluids())
            {
              case 1:
                {
                  const auto diffusivity_model_fluid =
                    properties_manager.get_tracer_diffusivity();
                  const auto diffusivity_model_solid =
                    properties_manager.get_tracer_diffusivity(0, 1);

                  diffusivity_model_fluid->vector_value(fields,
                                                        tracer_diffusivity);
                  diffusivity_model_solid->vector_value(fields,
                                                        tracer_diffusivity_1);

                  // We let the solid diffusivity model manage diffusivity
                  // assignation close to the solid
                  bool at_least_one_solid_q_point = false;
                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    if (sdf_values[q] < 0)
                      {
                        at_least_one_solid_q_point = true;
                        break;
                      }

                  if (at_least_one_solid_q_point)
                    for (unsigned int q = 0; q < this->n_q_points; ++q)
                      tracer_diffusivity[q] = tracer_diffusivity_1[q];
                  break;
                }
              default:
                throw std::runtime_error(
                  "Unsupported number of fluids and solids (1 solid and >1 fluid)");
            }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of solids (>1)");
    }
}


template class TracerScratchData<2>;
template class TracerScratchData<3>;
