#include <core/bdf.h>
#include <core/sdirk.h>

#include <solvers/navier_stokes_scratch_data.h>

#include <dem/dem.h>
#include <dem/dem_properties.h>

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

template <int dim>
void
NavierStokesScratchData<dim>::enable_free_surface(
  const FiniteElement<dim> &fe,
  const Quadrature<dim> &   quadrature,
  const Mapping<dim> &      mapping)
{
  gather_free_surface    = true;
  fe_values_free_surface = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Free surface
  phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_void_fraction(
  const FiniteElement<dim> &fe,
  const Quadrature<dim> &   quadrature,
  const Mapping<dim> &      mapping)
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
  const unsigned int n_global_max_particles_per_cell)
{
  gather_particles_information     = true;
  max_number_of_particles_per_cell = n_global_max_particles_per_cell;

  // Velocities
  particle_velocity =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  fluid_velocity_at_particle_location =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  cell_void_fraction = std::vector<double>(n_global_max_particles_per_cell);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_heat_transfer(
  const FiniteElement<dim> &fe,
  const Quadrature<dim> &   quadrature,
  const Mapping<dim> &      mapping)
{
  gather_temperature = true;
  fe_values_temperature =
    std::make_shared<FEValues<dim>>(mapping, fe, quadrature, update_values);

  temperature_values = std::vector<double>(this->n_q_points);
}


template class NavierStokesScratchData<2>;
template class NavierStokesScratchData<3>;
