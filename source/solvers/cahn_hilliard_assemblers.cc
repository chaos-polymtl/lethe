#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/cahn_hilliard_assemblers.h>
#include <solvers/copy_data.h>


template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Gather physical properties
  const double well_height       = this->ch_parameters.well_height;
  const double mobility_constant = this->ch_parameters.mobility_constant;
  const double epsilon           = scratch_data.epsilon;

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  auto &local_matrix = copy_data.local_matrix;

  // Constant mobility model
  if (this->ch_parameters.mobility_model == Parameters::MobilityModel::constant)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Gather into local variables the relevant fields
          const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

          // Store JxW in local variable for faster access;
          const double JxW = JxW_vec[q];

          const double phase_order_value = scratch_data.phase_order_values[q];

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
                     // Second equation
                     phi_potential_i * phi_potential_j -
                     4 * well_height * phi_potential_i *
                       (3 * phase_order_value * phase_order_value - 1.0) *
                       phi_phase_j -
                     epsilon * epsilon * grad_phi_potential_i *
                       grad_phi_phase_j) *
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
          const double phase_order_value = scratch_data.phase_order_values[q];

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
                     // Second equation
                     phi_potential_i * phi_potential_j -
                     4 * well_height * phi_potential_i *
                       (3 * phase_order_value * phase_order_value - 1.0) *
                       phi_phase_j -
                     epsilon * epsilon * grad_phi_potential_i *
                       grad_phi_phase_j) *
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
  StabilizedMethodsCopyData &   copy_data)
{
  // Gather physical properties
  const double well_height       = this->ch_parameters.well_height;
  const double mobility_constant = this->ch_parameters.mobility_constant;
  const double epsilon           = scratch_data.epsilon;

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  auto &local_rhs = copy_data.local_rhs;

  // Constant mobility model
  if (this->ch_parameters.mobility_model == Parameters::MobilityModel::constant)
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
          const double source_phase_order = scratch_data.source_phase_order[q];
          const double source_chemical_potential =
            scratch_data.source_chemical_potential[q];

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
                 // Second equation
                 - phi_potential_i * potential_value +
                 4 * well_height * phi_potential_i *
                   (phase_order_value * phase_order_value - 1) *
                   phase_order_value +
                 epsilon * epsilon * grad_phi_potential_i * phase_order_gradient
                 // Source term
                 + source_phase_order * phi_phase_i +
                 source_chemical_potential * phi_potential_i) *
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

          const double source_phase_order = scratch_data.source_phase_order[q];
          const double source_chemical_potential =
            scratch_data.source_chemical_potential[q];

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
                 // Second equation
                 - phi_potential_i * potential_value +
                 4 * well_height * phi_potential_i *
                   (phase_order_value * phase_order_value - 1) *
                   phase_order_value +
                 epsilon * epsilon * grad_phi_potential_i * phase_order_gradient
                 // Source term
                 + source_phase_order * phi_phase_i +
                 source_chemical_potential * phi_potential_i) *
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
  StabilizedMethodsCopyData &   copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_matrix = copy_data.local_matrix;

  for (unsigned int i_bc = 0; i_bc < this->boundary_conditions_ch.size; ++i_bc)
    {
      if (this->boundary_conditions_ch.type[i_bc] ==
          BoundaryConditions::BoundaryType::ch_angle_of_contact)
        {
          const double angle_of_contact =
            this->boundary_conditions_ch.angle_of_contact[i_bc];
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_ch.id[i_bc])
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
  StabilizedMethodsCopyData &   copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const double epsilon = scratch_data.epsilon;

  auto &local_rhs = copy_data.local_rhs;

  for (unsigned int i_bc = 0; i_bc < this->boundary_conditions_ch.size; ++i_bc)
    {
      if (this->boundary_conditions_ch.type[i_bc] ==
          BoundaryConditions::BoundaryType::ch_angle_of_contact)
        {
          const double angle_of_contact =
            this->boundary_conditions_ch.angle_of_contact[i_bc];
          for (unsigned int f = 0; f < scratch_data.n_faces; f++)
            {
              if (scratch_data.boundary_face_id[f] ==
                  this->boundary_conditions_ch.id[i_bc])
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
CahnHilliardAssemblerBDF<dim>::assemble_matrix(
  CahnHilliardScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &   copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

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
  StabilizedMethodsCopyData &   copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

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
