#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_assemblers.h>
#include <solvers/copy_data.h>



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Gather physical properties
  const double well_height = this->cahn_hilliard_parameters.well_height;
  const double mobility_constant =
    this->cahn_hilliard_parameters.cahn_hilliard_mobility_constant;
  const auto mobility_model =
    this->cahn_hilliard_parameters.cahn_hilliard_mobility_model;
  // std::cout<< "mobility = "<< mobility_constant<<std::endl;
  const double epsilon   = scratch_data.epsilon;
  const double cell_size = scratch_data.cell_size;
  const double xi =
    this->cahn_hilliard_parameters.potential_smoothing_coefficient;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Constant mobility model
  if (mobility_model == Parameters::CahnHilliardMobilityModel::constant)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Gather into local variables the relevant fields
          const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

          // Store JxW in local variable for faster access;
          const double JxW = JxW_vec[q];

          const double phase_order_value = scratch_data.phase_order_values[q];

          // Calculation of the magnitude of the velocity for the
          // stabilization parameter and the compression term for the phase
          // indicator
          const double u_mag = std::max(velocity_field.norm(), 1e-12);

          // Calculation of the GLS stabilization parameter. The
          // stabilization parameter used is different if the simulation is
          // steady or unsteady. In the unsteady case it includes the value
          // of the time-step. Hypothesis : advection dominated problem
          // (Pe>3) [Bochev et al., Stability of the SUPG finite element
          // method for transient advection-diffusion problems, CMAME 2004]
          const double tau =
            is_steady(method) ?
              1. / std::sqrt(std::pow(2. * u_mag / cell_size, 2) +
                             9 * std::pow(4 * mobility_constant /
                                            (cell_size * cell_size),
                                          2)) :
              1. /
                std::sqrt(
                  std::pow(sdt, 2) + std::pow(2. * u_mag / cell_size, 2) +
                  9 * std::pow(4 * mobility_constant / (cell_size * cell_size),
                               2));

          // Compute the strong jacobian vector
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> grad_phi_phase_j =
                scratch_data.grad_phi_potential[q][j];

              const double laplacian_phi_potential_j =
                scratch_data.laplacian_phi_potential[q][j];

              // Strong Jacobian associated with the GLS
              // stabilization
              strong_jacobian_vec[q][j] +=
                velocity_field * grad_phi_phase_j -
                mobility_constant * laplacian_phi_potential_j;
            }

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double         phi_phase_i = scratch_data.phi_phase[q][i];
              const Tensor<1, dim> grad_phi_phase_i =
                scratch_data.grad_phi_phase[q][i];

              const double phi_potential_i = scratch_data.phi_potential[q][i];
              const Tensor<1, dim> grad_phi_potential_i =
                scratch_data.grad_phi_potential[q][i];


              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const double phi_phase_j = scratch_data.phi_phase[q][j];
                  const Tensor<1, dim> grad_phi_phase_j =
                    scratch_data.grad_phi_phase[q][j];

                  const double phi_potential_j =
                    scratch_data.phi_potential[q][j];
                  const Tensor<1, dim> grad_phi_potential_j =
                    scratch_data.grad_phi_potential[q][j];

                  local_matrix(i, j) +=
                    // First equation
                    (phi_phase_i * (velocity_field * grad_phi_phase_j) +
                     mobility_constant * grad_phi_phase_i *
                       grad_phi_potential_j +
                     phi_potential_i * phi_potential_j -
                     // Second equation (Lovric et al.)
                     4 * well_height * phi_potential_i *
                       (3 * phase_order_value * phase_order_value - 1.0) *
                       phi_phase_j -
                     epsilon * epsilon * grad_phi_potential_i * grad_phi_phase_j
                     // Chemical potential smoothing
                     + xi * cell_size * cell_size * grad_phi_potential_i *
                         grad_phi_potential_j) *
                    JxW;

                  // Addition to the cell matrix for GLS stabilization
                  local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                        (grad_phi_phase_i * velocity_field) *
                                        JxW;
                }
            }
        }
    } // end loop on quadrature points
  // Quartic mobility model
  else
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Gather into local variables the relevant fields
          const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

          // Store JxW in local variable for faster access;
          const double         JxW = JxW_vec[q];
          const Tensor<1, dim> potential_gradient =
            scratch_data.chemical_potential_gradients[q];
          const double potential_laplacian =
            scratch_data.chemical_potential_laplacians[q];
          const double phase_order_value = scratch_data.phase_order_values[q];
          const Tensor<1, dim> phase_gradient =
            scratch_data.phase_order_gradients[q];

          // Calculation of the magnitude of the velocity for the
          // stabilization parameter and the compression term for the phase
          // indicator
          const double u_mag = std::max(velocity_field.norm(), 1e-12);

          // Calculation of the GLS stabilization parameter. The
          // stabilization parameter used is different if the simulation is
          // steady or unsteady. In the unsteady case it includes the value
          // of the time-step. Hypothesis : advection dominated problem
          // (Pe>3) [Bochev et al., Stability of the SUPG finite element
          // method for transient advection-diffusion problems, CMAME 2004]
          const double tau =
            is_steady(method) ?
              1. / std::sqrt(std::pow(2. * u_mag / cell_size, 2) +
                             9 * std::pow(4 * mobility_constant /
                                            (cell_size * cell_size),
                                          2)) :
              1. /
                std::sqrt(
                  std::pow(sdt, 2) + std::pow(2. * u_mag / cell_size, 2) +
                  9 * std::pow(4 * mobility_constant / (cell_size * cell_size),
                               2));

          // Compute the strong jacobian vector
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double         phi_phase_j = scratch_data.phi_phase[q][j];
              const Tensor<1, dim> grad_phi_phase_j =
                scratch_data.grad_phi_potential[q][j];
              const Tensor<1, dim> grad_phi_potential_j =
                scratch_data.grad_phi_potential[q][j];
              const double laplacian_phi_potential_j =
                scratch_data.laplacian_phi_potential[q][j];


              // Strong Jacobian associated with the GLS
              // stabilization
              strong_jacobian_vec[q][j] +=
                velocity_field * grad_phi_phase_j +
                4 * mobility_constant * phi_phase_j *
                  ((1 - phase_order_value * phase_order_value) *
                     phase_gradient * potential_gradient -
                   2 * phase_order_value * phase_order_value * phase_gradient *
                     potential_gradient +
                   phase_order_value *
                     (1 - phase_order_value * phase_order_value) *
                     potential_laplacian) +
                4 * mobility_constant * phase_order_value *
                  (1 - phase_order_value * phase_order_value) *
                  grad_phi_phase_j * potential_gradient +
                4 * mobility_constant * phase_order_value *
                  (1 - phase_order_value * phase_order_value) * phase_gradient *
                  grad_phi_potential_j -
                mobility_constant *
                  (1 - phase_order_value * phase_order_value) *
                  (1 - phase_order_value * phase_order_value) *
                  laplacian_phi_potential_j;
            }


          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double         phi_phase_i = scratch_data.phi_phase[q][i];
              const Tensor<1, dim> grad_phi_phase_i =
                scratch_data.grad_phi_phase[q][i];

              const double phi_potential_i = scratch_data.phi_potential[q][i];
              const Tensor<1, dim> grad_phi_potential_i =
                scratch_data.grad_phi_potential[q][i];

              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const double phi_phase_j = scratch_data.phi_phase[q][j];
                  const Tensor<1, dim> grad_phi_phase_j =
                    scratch_data.grad_phi_phase[q][j];

                  const double phi_potential_j =
                    scratch_data.phi_potential[q][j];
                  const Tensor<1, dim> grad_phi_potential_j =
                    scratch_data.grad_phi_potential[q][j];

                  local_matrix(i, j) +=
                    // First equation
                    (phi_phase_i * (velocity_field * grad_phi_phase_j) -
                     4 * mobility_constant * grad_phi_phase_i *
                       potential_gradient * phase_order_value *
                       (1 - phase_order_value * phase_order_value) *
                       phi_phase_j +
                     mobility_constant * grad_phi_phase_i *
                       grad_phi_potential_j *
                       (1 - phase_order_value * phase_order_value) *
                       (1 - phase_order_value * phase_order_value) +
                     phi_potential_i * phi_potential_j -
                     // Second equation (Lovric et al.)
                     4 * well_height * phi_potential_i *
                       (3 * phase_order_value * phase_order_value - 1.0) *
                       phi_phase_j -
                     epsilon * epsilon * grad_phi_potential_i * grad_phi_phase_j
                     // Chemical potential smoothing
                     + xi * cell_size * cell_size * grad_phi_potential_i *
                         grad_phi_potential_j) *
                    JxW;

                  // Addition to the cell matrix for GLS stabilization
                  local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                        (grad_phi_phase_i * velocity_field) *
                                        JxW;
                }
            }
        }
    } // end loop on quadrature points
}



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_rhs(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Gather physical properties
  const double well_height = this->cahn_hilliard_parameters.well_height;
  const double mobility_constant =
    this->cahn_hilliard_parameters.cahn_hilliard_mobility_constant;
  const auto mobility_model =
    this->cahn_hilliard_parameters.cahn_hilliard_mobility_model;
  const double epsilon   = scratch_data.epsilon;
  const double cell_size = scratch_data.cell_size;
  const double xi =
    this->cahn_hilliard_parameters.potential_smoothing_coefficient;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Constant mobility model
  if (mobility_model == Parameters::CahnHilliardMobilityModel::constant)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Gather into local variables the relevant fields
          const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

          // Store JxW in local variable for faster access;
          const double JxW               = JxW_vec[q];
          const double phase_order_value = scratch_data.phase_order_values[q];
          const Tensor<1, dim> phase_order_gradient =
            scratch_data.phase_order_gradients[q];
          const double potential_value =
            scratch_data.chemical_potential_values[q];
          const Tensor<1, dim> potential_gradient =
            scratch_data.chemical_potential_gradients[q];
          const double potential_laplacian =
            scratch_data.chemical_potential_laplacians[q];
          const double source_phase_order = scratch_data.source_phase_order[q];
          const double source_chemical_potential =
            scratch_data.source_chemical_potential[q];

          // Calculation of the magnitude of the velocity for the
          // stabilization parameter
          const double u_mag = std::max(velocity_field.norm(), 1e-12);

          // Calculation of the GLS stabilization parameter. The
          // stabilization parameter used is different if the simulation is
          // steady or unsteady. In the unsteady case it includes the value
          // of the time-step. Hypothesis : advection dominated problem
          // (Pe>3) [Bochev et al., Stability of the SUPG finite element
          // method for transient advection-diffusion problems, CMAME 2004]
          const double tau =
            is_steady(method) ?
              1. / std::sqrt(std::pow(2. * u_mag / cell_size, 2) +
                             9 * std::pow(4 * mobility_constant /
                                            (cell_size * cell_size),
                                          2)) :
              1. /
                std::sqrt(
                  std::pow(sdt, 2) + std::pow(2. * u_mag / cell_size, 2) +
                  9 * std::pow(4 * mobility_constant / (cell_size * cell_size),
                               2));

          // Calculate the strong residual for GLS stabilization
          strong_residual_vec[q] += velocity_field * phase_order_gradient -
                                    mobility_constant * potential_laplacian;

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double         phi_phase_i = scratch_data.phi_phase[q][i];
              const Tensor<1, dim> grad_phi_phase_i =
                scratch_data.grad_phi_phase[q][i];
              const double phi_potential_i = scratch_data.phi_potential[q][i];
              const Tensor<1, dim> grad_phi_potential_i =
                scratch_data.grad_phi_potential[q][i];

              local_rhs(i) +=
                // First equation
                (-phi_phase_i * (velocity_field * phase_order_gradient) -
                 mobility_constant * grad_phi_phase_i * potential_gradient
                 // Second equation (2nd article)
                 - phi_potential_i * potential_value +
                 // Second equation (Lovric et al.)
                 4 * well_height * phi_potential_i *
                   (phase_order_value * phase_order_value - 1) *
                   phase_order_value +
                 epsilon * epsilon * grad_phi_potential_i * phase_order_gradient
                 // Chemical potential smoothing
                 - xi * cell_size * cell_size * grad_phi_potential_i *
                     potential_gradient
                 // Source term
                 + source_phase_order * phi_phase_i +
                 source_chemical_potential * phi_potential_i) *
                JxW;

              // Addition to the RHS for GLS stabilization
              local_rhs(i) +=
                -tau *
                (strong_residual_vec[q] * (grad_phi_phase_i * velocity_field)) *
                JxW;
            }
        }
    } // end loop on quadrature points
  // Quartic mobility model
  else
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Gather into local variables the relevant fields
          const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

          // Store JxW in local variable for faster access;
          const double JxW               = JxW_vec[q];
          const double phase_order_value = scratch_data.phase_order_values[q];
          const Tensor<1, dim> phase_order_gradient =
            scratch_data.phase_order_gradients[q];
          const double potential_value =
            scratch_data.chemical_potential_values[q];
          const Tensor<1, dim> potential_gradient =
            scratch_data.chemical_potential_gradients[q];
          const double potential_laplacian =
            scratch_data.chemical_potential_laplacians[q];

          const double source_phase_order = scratch_data.source_phase_order[q];
          const double source_chemical_potential =
            scratch_data.source_chemical_potential[q];

          // Calculation of the magnitude of the velocity for the
          // stabilization parameter
          const double u_mag = std::max(velocity_field.norm(), 1e-12);

          // Calculation of the GLS stabilization parameter. The
          // stabilization parameter used is different if the simulation is
          // steady or unsteady. In the unsteady case it includes the value
          // of the time-step. Hypothesis : advection dominated problem
          // (Pe>3) [Bochev et al., Stability of the SUPG finite element
          // method for transient advection-diffusion problems, CMAME 2004]
          const double tau =
            is_steady(method) ?
              1. / std::sqrt(std::pow(2. * u_mag / cell_size, 2) +
                             9 * std::pow(4 * mobility_constant /
                                            (cell_size * cell_size),
                                          2)) :
              1. /
                std::sqrt(
                  std::pow(sdt, 2) + std::pow(2. * u_mag / cell_size, 2) +
                  9 * std::pow(4 * mobility_constant / (cell_size * cell_size),
                               2));

          // Calculate the strong residual for GLS stabilization
          strong_residual_vec[q] +=
            velocity_field * phase_order_gradient +
            4 * mobility_constant * phase_order_value *
              (1 - phase_order_value * phase_order_value) *
              phase_order_gradient * potential_gradient -
            mobility_constant * (1 - phase_order_value * phase_order_value) *
              (1 - phase_order_value * phase_order_value) * potential_laplacian;

          for (unsigned int i = 0; i < n_dofs; ++i)
            {
              const double         phi_phase_i = scratch_data.phi_phase[q][i];
              const Tensor<1, dim> grad_phi_phase_i =
                scratch_data.grad_phi_phase[q][i];
              const double phi_potential_i = scratch_data.phi_potential[q][i];
              const Tensor<1, dim> grad_phi_potential_i =
                scratch_data.grad_phi_potential[q][i];

              local_rhs(i) +=
                // First equation
                (-phi_phase_i * (velocity_field * phase_order_gradient) -
                 grad_phi_phase_i * potential_gradient * mobility_constant *
                   (1 - phase_order_value * phase_order_value) *
                   (1 - phase_order_value * phase_order_value)
                 // Second equation (2nd article)
                 - phi_potential_i * potential_value +
                 // Second equation (Lovric et al.)
                 4 * well_height * phi_potential_i *
                   (phase_order_value * phase_order_value - 1) *
                   phase_order_value +
                 epsilon * epsilon * grad_phi_potential_i * phase_order_gradient
                 // Chemical potential smoothing
                 - xi * cell_size * cell_size * grad_phi_potential_i *
                     potential_gradient
                 // Source term
                 + source_phase_order * phi_phase_i +
                 source_chemical_potential * phi_potential_i) *
                JxW;

              // Addition to the RHS for GLS stabilization
              local_rhs(i) +=
                -tau *
                (strong_residual_vec[q] * (grad_phi_phase_i * velocity_field)) *
                JxW;
            }
        }
    } // end loop on quadrature points
}

template class CahnHilliardAssemblerCore<2>;
template class CahnHilliardAssemblerCore<3>;

template <int dim>
void
CahnHilliardAssemblerAngleOfContact<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_matrix = copy_data.local_matrix;

  for (unsigned int i_bc = 0;
       i_bc < this->boundary_conditions_cahn_hilliard.size;
       ++i_bc)
    {
      if (this->boundary_conditions_cahn_hilliard.type[i_bc] ==
          BoundaryConditions::BoundaryType::cahn_hilliard_angle_of_contact)
        {
          const double angle_of_contact =
            this->boundary_conditions_cahn_hilliard.angle_of_contact[i_bc];
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_cahn_hilliard.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const Tensor<1, dim> face_phase_grad_value =
                        scratch_data.face_phase_grad_values[f][q];

                      const double JxW_face = scratch_data.face_JxW[f][q];

                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_potential_i =
                            scratch_data.phi_potential[q][i];

                          for (unsigned int j = 0; j < scratch_data.n_dofs; ++j)
                            {
                              const Tensor<1, dim> grad_phi_face_phase_j =
                                scratch_data.grad_phi_face_phase[f][q][j];

                              // 1st article : Lovric et al.
                              local_matrix(i, j) +=
                                -epsilon * epsilon * phi_potential_i *
                                grad_phi_face_phase_j * face_phase_grad_value *
                                std::cos(angle_of_contact * M_PI / 180.0) *
                                (1.0 / (face_phase_grad_value.norm() +
                                        std::numeric_limits<double>::min())) *
                                JxW_face;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
CahnHilliardAssemblerAngleOfContact<dim>::assemble_rhs(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_rhs = copy_data.local_rhs;

  for (unsigned int i_bc = 0;
       i_bc < this->boundary_conditions_cahn_hilliard.size;
       ++i_bc)
    {
      if (this->boundary_conditions_cahn_hilliard.type[i_bc] ==
          BoundaryConditions::BoundaryType::cahn_hilliard_angle_of_contact)
        {
          const double angle_of_contact =
            this->boundary_conditions_cahn_hilliard.angle_of_contact[i_bc];
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_cahn_hilliard.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const Tensor<1, dim> face_phase_grad_value =
                        scratch_data.face_phase_grad_values[f][q];

                      const double JxW_face = scratch_data.face_JxW[f][q];

                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_potential_i =
                            scratch_data.phi_potential[q][i];

                          local_rhs(i) +=
                            epsilon * epsilon * phi_potential_i *
                            std::cos(angle_of_contact * M_PI / 180.0) *
                            (face_phase_grad_value.norm()) * JxW_face;
                        }
                    }
                }
            }
        }
    }
}

template class CahnHilliardAssemblerAngleOfContact<2>;
template class CahnHilliardAssemblerAngleOfContact<3>;

template <int dim>
void
CahnHilliardAssemblerFreeAngle<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_matrix = copy_data.local_matrix;

  for (unsigned int i_bc = 0;
       i_bc < this->boundary_conditions_cahn_hilliard.size;
       ++i_bc)
    {
      if (this->boundary_conditions_cahn_hilliard.type[i_bc] ==
          BoundaryConditions::BoundaryType::cahn_hilliard_free_angle)
        {
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_cahn_hilliard.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const Tensor<1, dim> face_normal =
                        scratch_data.face_normal[f][q];
                      const double JxW_face = scratch_data.face_JxW[f][q];

                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_potential_i =
                            scratch_data.phi_potential[q][i];

                          for (unsigned int j = 0; j < scratch_data.n_dofs; ++j)
                            {
                              const Tensor<1, dim> grad_phi_face_phase_j =
                                scratch_data.grad_phi_face_phase[f][q][j];
                              // 2nd article : Lovric et al.
                              local_matrix(i, j) -=
                                epsilon * epsilon * phi_potential_i *
                                grad_phi_face_phase_j * face_normal * JxW_face;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void
CahnHilliardAssemblerFreeAngle<dim>::assemble_rhs(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_rhs = copy_data.local_rhs;

  for (unsigned int i_bc = 0;
       i_bc < this->boundary_conditions_cahn_hilliard.size;
       ++i_bc)
    {
      if (this->boundary_conditions_cahn_hilliard.type[i_bc] ==
          BoundaryConditions::BoundaryType::cahn_hilliard_free_angle)
        {
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_cahn_hilliard.id[i_bc])
                {
                  for (unsigned int q = 0; q < scratch_data.n_faces_q_points;
                       ++q)
                    {
                      const Tensor<1, dim> face_phase_grad_value =
                        scratch_data.face_phase_grad_values[f][q];
                      const Tensor<1, dim> face_normal =
                        scratch_data.face_normal[f][q];
                      const double JxW_face = scratch_data.face_JxW[f][q];

                      for (unsigned int i = 0; i < scratch_data.n_dofs; ++i)
                        {
                          const double phi_potential_i =
                            scratch_data.phi_potential[q][i];

                          local_rhs(i) += epsilon * epsilon * phi_potential_i *
                                          face_phase_grad_value * face_normal *
                                          JxW_face;
                        }
                    }
                }
            }
        }
    }
}

template class CahnHilliardAssemblerFreeAngle<2>;
template class CahnHilliardAssemblerFreeAngle<3>;


template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &strong_residual = copy_data.strong_residual;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> phase_order(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_order[0] = scratch_data.phase_order_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_order[p + 1] = scratch_data.previous_phase_order_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_order[p]);
        }
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi_phase[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_phase_i = scratch_data.phi_phase[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double phi_phase_j = scratch_data.phi_phase[q][j];

              local_matrix(i, j) +=
                phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_rhs(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto          method = this->simulation_control->get_assembly_method();
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();

  // Vector for the BDF coefficients
  Vector<double>      bdf_coefs = bdf_coefficients(method, time_steps_vector);
  std::vector<double> phase_order(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_order[0] = scratch_data.phase_order_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_order[p + 1] = scratch_data.previous_phase_order_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_order[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_phase_i = scratch_data.phi_phase[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (phase_order[p] * phi_phase_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class CahnHilliardAssemblerBDF<2>;
template class CahnHilliardAssemblerBDF<3>;
