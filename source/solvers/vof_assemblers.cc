// SPDX-FileCopyrightText: Copyright (c) 2021-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/vof_assemblers.h>

#include <deal.II/fe/fe_interface_values.h>

template <int dim>
void
VOFAssemblerCore<dim>::assemble_matrix(const VOFScratchData<dim> &scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Add a small diffusivity, used to artificially diffuse VOF
  const double diffusivity = this->vof_parameters.diffusivity;

  // Copy data elements
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Assemble local matrix and right-hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity    = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];

      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      // Define tolerance to avoid division by zero
      const double tolerance = 1e-12;

      // Compute the velocity magnitude for the stabilization parameter and the
      // compression term for the phase indicator
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case, it includes the value
      // of the time-step. Hypothesis: advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(
                 Utilities::fixed_power<2>(2. * u_mag / h) +
                 9 * Utilities::fixed_power<2>(4 * diffusivity / (h * h))) :
          1. /
            std::sqrt(Utilities::fixed_power<2>(sdt) +
                      Utilities::fixed_power<2>(2. * u_mag / h) +
                      9 * Utilities::fixed_power<2>(4 * diffusivity / (h * h)));

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto phi_phase_j           = scratch_data.phi[q][j];
          const auto grad_phi_phase_j      = scratch_data.grad_phi[q][j];
          const auto laplacian_phi_phase_j = scratch_data.laplacian_phi[q][j];

          // Strong Jacobian associated with the GLS stabilization
          strong_jacobian_vec[q][j] +=
            velocity * grad_phi_phase_j - diffusivity * laplacian_phi_phase_j;

          if (compressible)
            strong_jacobian_vec[q][j] += phi_phase_j * velocity_divergence;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j      = scratch_data.phi[q][j];
              const auto grad_phi_phase_j = scratch_data.grad_phi[q][j];

              // Weak form for advection-diffusion:
              // u * grad(phase) + phase * grad(u) - diffusivity *
              // laplacian(phase) = 0
              local_matrix(i, j) +=
                (phi_phase_i * velocity * grad_phi_phase_j +
                 diffusivity * grad_phi_phase_i * grad_phi_phase_j) *
                JxW;

              // Add compressibility term if the VOF is compressible
              if (compressible)
                local_matrix(i, j) +=
                  phi_phase_i * phi_phase_j * velocity_divergence * JxW;

              // Addition to the cell matrix for GLS stabilization
              local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                    (grad_phi_phase_i * velocity) * JxW;
            }
        }
    } // end loop on quadrature points
}


template <int dim>
void
VOFAssemblerCore<dim>::assemble_rhs(const VOFScratchData<dim> &scratch_data,
                                    StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto method = this->simulation_control->get_assembly_method();

  // Add a small diffusivity, used to artificially diffuse VOF
  const double diffusivity = this->vof_parameters.diffusivity;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Time steps and inverse time steps which is used for stabilization constant
  std::vector<double> time_steps_vector =
    this->simulation_control->get_time_steps_vector();
  const double dt  = time_steps_vector[0];
  const double sdt = 1. / dt;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Assemble local matrix and right-hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         phase = scratch_data.present_phase_values[q];
      const Tensor<1, dim> phase_gradient   = scratch_data.phase_gradients[q];
      const double         phase_laplacians = scratch_data.phase_laplacians[q];
      const Tensor<1, dim> velocity         = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Define tolerance to avoid division by zero
      const double tolerance = 1e-12;

      // Compute the velocity magnitude for the stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case, it includes the value
      // of the time-step. Hypothesis: advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(
                 Utilities::fixed_power<2>(2. * u_mag / h) +
                 9 * Utilities::fixed_power<2>(4 * diffusivity / (h * h))) :
          1. /
            std::sqrt(Utilities::fixed_power<2>(sdt) +
                      Utilities::fixed_power<2>(2. * u_mag / h) +
                      9 * Utilities::fixed_power<2>(4 * diffusivity / (h * h)));

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] +=
        velocity * phase_gradient - diffusivity * phase_laplacians;

      if (compressible)
        strong_residual_vec[q] += phase * velocity_divergence;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];

          // rhs for: u * grad(phase) + phase * grad(u) - diffusivity *
          // laplacian(phase) = 0
          local_rhs(i) -= (phi_phase_i * velocity * phase_gradient +
                           diffusivity * grad_phi_phase_i * phase_gradient) *
                          JxW;

          if (compressible)
            local_rhs(i) -= phi_phase_i * phase * velocity_divergence * JxW;

          local_rhs(i) -=
            tau * (strong_residual_vec[q] * (grad_phi_phase_i * velocity)) *
            JxW;
        }
    }
}

template class VOFAssemblerCore<2>;
template class VOFAssemblerCore<3>;


template <int dim>
void
VOFAssemblerBDF<dim>::assemble_matrix(const VOFScratchData<dim> &scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
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
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_values[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] = scratch_data.previous_phase_values[p][q];


      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i = scratch_data.phi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j = scratch_data.phi[q][j];

              local_matrix(i, j) +=
                phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
VOFAssemblerBDF<dim>::assemble_rhs(const VOFScratchData<dim> &scratch_data,
                                   StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &local_rhs       = copy_data.local_rhs;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
  std::vector<double> phase_value(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      phase_value[0] = scratch_data.present_phase_values[q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        phase_value[p + 1] = scratch_data.previous_phase_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * phase_value[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_phase_i = scratch_data.phi[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (phase_value[p] * phi_phase_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class VOFAssemblerBDF<2>;
template class VOFAssemblerBDF<3>;


template <int dim>
void
VOFAssemblerDCDDStabilization<dim>::assemble_matrix(
  const VOFScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Assemble local matrix
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields for faster access
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];
      const double         JxW      = JxW_vec[q];

      // We use the previous phase gradient for the shock capture (explicit)
      const Tensor<1, dim> previous_phase_gradient =
        scratch_data.previous_phase_gradients[q];

      // Define a tolerance to avoid division by zero
      const double tolerance = 1e-12;

      // Compute the velocity magnitude and phase gradient norm for the
      // stabilization term
      const double velocity_norm       = velocity.norm();
      const double phase_gradient_norm = previous_phase_gradient.norm();

      // In Tezduyar 2003, this is denoted r in equation (67).
      Tensor<1, dim> gradient_unit_vector =
        previous_phase_gradient / (phase_gradient_norm + tolerance);

      // We remove the diffusion aligned with the velocity as is done in the
      // original article. In Tezduyar 2003, this is denoted s in equation (67).
      Tensor<1, dim> velocity_unit_vector =
        velocity / (velocity_norm + tolerance);

      // Directionality tensor describing [rr - (r⋅s)^2*ss]
      const Tensor<2, dim> dir_tensor =
        outer_product(gradient_unit_vector, gradient_unit_vector) -
        ((gradient_unit_vector * velocity_unit_vector) *
         outer_product(velocity_unit_vector, velocity_unit_vector));

      // Compute the artificial viscosity of the shock capture
      const double nu_dcdd =
        (0.5 * h * h) * velocity.norm() * phase_gradient_norm;

      // Assemble DCDD shock capturing contributions
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto grad_phi_phase_j = scratch_data.grad_phi[q][j];
              local_matrix(i, j) +=
                (nu_dcdd * scalar_product(grad_phi_phase_i,
                                          dir_tensor * grad_phi_phase_j)) *
                JxW;
            }
        }
    }
}

template <int dim>
void
VOFAssemblerDCDDStabilization<dim>::assemble_rhs(
  const VOFScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;
  const double       h          = scratch_data.cell_size;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Assemble local right-hand side (rhs)
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields for faster access
      const Tensor<1, dim> phase_gradient = scratch_data.phase_gradients[q];
      const Tensor<1, dim> velocity       = scratch_data.velocity_values[q];
      const double         JxW            = JxW_vec[q];

      // We use the previous phase gradient for the shock capture (explicit)
      const Tensor<1, dim> previous_phase_gradient =
        scratch_data.previous_phase_gradients[q];

      // Define a tolerance to avoid division by zero
      const double tolerance = 1e-12;

      // Compute the velocity magnitude and phase gradient norm for the
      // stabilization term
      const double velocity_norm       = velocity.norm();
      const double phase_gradient_norm = previous_phase_gradient.norm();

      // In Tezduyar 2003, this is denoted r in equation (67).
      Tensor<1, dim> gradient_unit_vector =
        previous_phase_gradient / (phase_gradient_norm + tolerance);

      // We remove the diffusion aligned with the velocity as is done in the
      // original article. In Tezduyar 2003, this is denoted s in equation (67).
      Tensor<1, dim> velocity_unit_vector =
        velocity / (velocity_norm + tolerance);

      // Directionality tensor describing [rr - (r⋅s)^2*ss]
      const Tensor<2, dim> dir_tensor =
        outer_product(gradient_unit_vector, gradient_unit_vector) -
        ((gradient_unit_vector * velocity_unit_vector) *
         outer_product(velocity_unit_vector, velocity_unit_vector));

      // Compute the artificial viscosity of the shock capture
      const double nu_dcdd =
        (0.5 * h * h) * velocity.norm() * phase_gradient_norm;

      // Assemble DCDD shock capturing contributions
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];
          local_rhs(i) -=
            nu_dcdd *
            scalar_product(grad_phi_phase_i, dir_tensor * phase_gradient) * JxW;
        }
    }
}

template class VOFAssemblerDCDDStabilization<2>;
template class VOFAssemblerDCDDStabilization<3>;


template <int dim>
void
VOFAssemblerInletOutlet<dim>::assemble_matrix(
  const VOFScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const FiniteElement<dim> &fe = scratch_data.fe_face_values_vof.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_matrix = copy_data.local_matrix;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];
          // Check if the face is part of a boundary on which we apply an
          // inlet-outlet
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::vof_inlet_outlet)
            {
              const double beta =
                this->boundary_conditions.beta.at(boundary_id);

              // Assemble the matrix of the BC
              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values_vof.dof_indices())
                    {
                      double normal_outflow =
                        (scratch_data.boundary_face_velocity_values[f][q] *
                         scratch_data.face_normal[f][q]);

                      for (const unsigned int j :
                           scratch_data.fe_face_values_vof.dof_indices())
                        {
                          if (normal_outflow < 0)
                            {
                              local_matrix(i, j) +=
                                -penalty_parameter * beta *
                                scratch_data.face_phi[f][q][j] *
                                scratch_data.face_phi[f][q][i] * JxW;
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
VOFAssemblerInletOutlet<dim>::assemble_rhs(
  const VOFScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &copy_data)
{
  if (!scratch_data.is_boundary_cell)
    return;

  const FiniteElement<dim> &fe = scratch_data.fe_face_values_vof.get_fe();

  const double penalty_parameter =
    1. / std::pow(scratch_data.cell_size, fe.degree + 1);
  auto &local_rhs = copy_data.local_rhs;

  // Loop over the faces of the cell.
  for (unsigned int f = 0; f < scratch_data.n_faces; ++f)
    {
      // Check if the face is on a boundary
      if (scratch_data.is_boundary_face[f])
        {
          types::boundary_id boundary_id = scratch_data.boundary_face_id[f];

          // Check if the face is part of the boundary that has a
          // weakly imposed Dirichlet BC.
          if (this->boundary_conditions.type.at(boundary_id) ==
              BoundaryConditions::BoundaryType::vof_inlet_outlet)
            {
              const double beta =
                this->boundary_conditions.beta.at(boundary_id);
              const double inlet_phase =
                this->boundary_conditions.inlet_phase.at(boundary_id);

              for (unsigned int q = 0; q < scratch_data.n_faces_q_points; ++q)
                {
                  const double JxW = scratch_data.face_JxW[f][q];
                  for (const unsigned int i :
                       scratch_data.fe_face_values_vof.dof_indices())
                    {
                      // Calculate beta term depending on the
                      // value of  u*n. If it is positive (outgoing
                      // flow) then no penalization is applied.
                      // If it is negative (inflow) then we penalize then
                      // we penalize the phase value to adhere to the
                      // prescribed inlet phase
                      const double normal_outflow =
                        (scratch_data.boundary_face_velocity_values[f][q] *
                         scratch_data.face_normal[f][q]);

                      if (normal_outflow < 0.)
                        {
                          local_rhs(i) +=
                            penalty_parameter * beta * normal_outflow *
                            ((scratch_data.face_phase_values[f][q] -
                              inlet_phase) *
                             scratch_data.face_phi[f][q][i]) *
                            JxW;
                        }
                    }
                }
            }
        }
    }
}

template class VOFAssemblerInletOutlet<2>;
template class VOFAssemblerInletOutlet<3>;



template <int dim>
void
VOFAssemblerDGCore<dim>::assemble_matrix(
  const VOFScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Assembling local matrix
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j = scratch_data.phi[q][j];

              // Weak form of : u * grad(phi) =0
              // Note that the advection term has been weakened for it to appear
              // explicitly in the weak form as a boundary term.
              local_matrix(i, j) +=
                (-grad_phi_phase_i * velocity * phi_phase_j) * JxW;
            }
        }
    } // end loop on quadrature points
}


template <int dim>
void
VOFAssemblerDGCore<dim>::assemble_rhs(const VOFScratchData<dim> &scratch_data,
                                      StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Assembling local right-hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         phase_value = scratch_data.present_phase_values[q];
      const Tensor<1, dim> velocity    = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];

          // Linearized weak form of the strong problem u * grad(phi) =0
          // Note that the advection term has been weakened for it to appear
          // explicitly in the weak form as a boundary term.
          local_rhs(i) -= (-grad_phi_phase_i * velocity * phase_value) * JxW;
        }
    } // end loop on quadrature points
}

template class VOFAssemblerDGCore<2>;
template class VOFAssemblerDGCore<3>;



template <int dim>
void
VOFAssemblerSIPG<dim>::assemble_matrix(const VOFScratchData<dim> &scratch_data,
                                       StabilizedDGMethodsCopyData &copy_data)
{
  const FEInterfaceValues<dim> &fe_iv    = scratch_data.fe_interface_values_vof;
  const auto                   &q_points = fe_iv.get_quadrature_points();
  auto                         &copy_data_face = copy_data.face_data.back();
  const unsigned int            n_dofs = fe_iv.n_current_interface_dofs();


  const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

  for (unsigned int q = 0; q < q_points.size(); ++q)
    {
      const double velocity_dot_n =
        scratch_data.face_velocity_values[q] * normals[q];
      for (unsigned int i = 0; i < n_dofs; ++i)
        for (unsigned int j = 0; j < n_dofs; ++j)
          {
            // Assembles symetric interior penalty method
            // There is no penalty coefficient here since the VOF equation
            // does not have diffusivity. The only flux that remain
            // is the upwinded advective flux.
            // ( [[φ]] ,u·n φ_upwind )
            copy_data_face.face_matrix(i, j) +=
              fe_iv.jump_in_shape_values(i, q) *
              fe_iv.shape_value((velocity_dot_n > 0.), j, q) * velocity_dot_n *
              JxW[q];
          }
    }
}

template <int dim>
void
VOFAssemblerSIPG<dim>::assemble_rhs(const VOFScratchData<dim>   &scratch_data,
                                    StabilizedDGMethodsCopyData &copy_data)
{
  const FEInterfaceValues<dim> &fe_iv    = scratch_data.fe_interface_values_vof;
  const auto                   &q_points = fe_iv.get_quadrature_points();
  const unsigned int            n_dofs   = fe_iv.n_current_interface_dofs();

  auto &copy_data_face = copy_data.face_data.back();

  const std::vector<double>         &JxW     = fe_iv.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();

  for (unsigned int q = 0; q < q_points.size(); ++q)
    {
      const double velocity_dot_n =
        scratch_data.face_velocity_values[q] * normals[q];
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          // Assemble advection terms with upwinding
          // If velocity_dot_n (u·n) > 0, then the upwind term needs to come
          // from here. Otherwise, it needs to come from there.
          if (velocity_dot_n > 0)
            {
              copy_data_face.face_rhs(i) -=
                fe_iv.jump_in_shape_values(i, q) // [\phi_i]
                * scratch_data.values_here[q]    // \phi_i^{upwind}
                * velocity_dot_n                 // (u . n)
                * JxW[q];                        // dx
            }
          else
            {
              copy_data_face.face_rhs(i) -=
                fe_iv.jump_in_shape_values(i, q) // [\phi_i]
                * scratch_data.values_there[q]   // \phi_i^{upwind}
                * velocity_dot_n                 // (u . n)
                * JxW[q];                        // dx
            }
        }
    }
}

template class VOFAssemblerSIPG<2>;
template class VOFAssemblerSIPG<3>;
