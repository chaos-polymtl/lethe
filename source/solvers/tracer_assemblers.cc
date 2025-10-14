// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/tracer_assemblers.h>


template <int dim>
void
TracerAssemblerCore<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &diffusivity_vector =
    scratch_data.tracer_diffusivity;
  const auto method = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
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
  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         diffusivity     = diffusivity_vector[q];
      const Tensor<1, dim> tracer_gradient = scratch_data.tracer_gradients[q];
      const Tensor<1, dim> velocity        = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // We use the previous concentration for the shock capture if the
      // simulation is transient
      const Tensor<1, dim> tracer_gradient_for_dcdd =
        is_steady(method) ? tracer_gradient :
                            scratch_data.previous_tracer_gradients[q];
      const double tracer_gradient_for_dcdd_norm =
        tracer_gradient_for_dcdd.norm();

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_tracer.get_fe().degree;

      const double vdcdd = (0.5 * h) * (velocity.norm() * velocity.norm()) *
                           pow(tracer_gradient_for_dcdd.norm() * h, order);

      const double   tolerance = 1e-12;
      Tensor<1, dim> s         = velocity / (velocity.norm() + tolerance);
      Tensor<1, dim> r =
        tracer_gradient_for_dcdd / (tracer_gradient_for_dcdd_norm + tolerance);

      const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> rr          = outer_product(r, r);
      const Tensor<2, dim> dcdd_factor = rr - k_corr;

      // If the method is not steady, we use the previous gradient to calculate
      // the DCDD schock capture term and consequently we do not need an
      // estimation of the derivative of the vdcdd term.
      const double d_vdcdd =
        is_steady(method) ?
          order * (0.5 * h * h) * (velocity.norm() * velocity.norm()) *
            pow(tracer_gradient_for_dcdd_norm * h, order - 1) :
          0;



      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi[q][j];
          const double laplacian_phi_T_j    = scratch_data.laplacian_phi[q][j];
          strong_jacobian_vec[q][j] +=
            velocity * grad_phi_T_j - diffusivity * laplacian_phi_T_j;
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi[q][j];

              // Weak form : - D * laplacian T +  u * gradT - f=0
              local_matrix(i, j) += (diffusivity * grad_phi_T_i * grad_phi_T_j +
                                     phi_T_i * velocity * grad_phi_T_j) *
                                    JxW;



              local_matrix(i, j) += tau * strong_jacobian_vec[q][j] *
                                    (grad_phi_T_i * velocity) * JxW;

              local_matrix(i, j) +=
                (vdcdd *
                   scalar_product(grad_phi_T_i, dcdd_factor * grad_phi_T_j) +
                 d_vdcdd * grad_phi_T_j.norm() *
                   scalar_product(grad_phi_T_i,
                                  dcdd_factor * tracer_gradient)) *
                JxW;
            }
        }
    } // end loop on quadrature points
}


template <int dim>
void
TracerAssemblerCore<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &diffusivity_vector =
    scratch_data.tracer_diffusivity;
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

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         diffusivity     = diffusivity_vector[q];
      const Tensor<1, dim> tracer_gradient = scratch_data.tracer_gradients[q];
      const Tensor<1, dim> previous_tracer_gradient =
        scratch_data.previous_tracer_gradients[q];
      const double         tracer_laplacian = scratch_data.tracer_laplacians[q];
      const Tensor<1, dim> velocity         = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // We use the previous concentration for the shock capture if the
      // simulation is transient
      const Tensor<1, dim> tracer_gradient_for_dcdd =
        is_steady(method) ? tracer_gradient : previous_tracer_gradient;
      const double tracer_gradient_for_dcdd_norm =
        tracer_gradient_for_dcdd.norm();

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_tracer.get_fe().degree;


      const double vdcdd = (0.5 * h) * (velocity.norm() * velocity.norm()) *
                           pow(tracer_gradient_for_dcdd_norm * h, order);

      const double   tolerance = 1e-12;
      Tensor<1, dim> s         = velocity / (velocity.norm() + tolerance);
      Tensor<1, dim> r =
        tracer_gradient_for_dcdd / (tracer_gradient_for_dcdd_norm + tolerance);

      const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> rr          = outer_product(r, r);
      const Tensor<2, dim> dcdd_factor = rr - k_corr;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step
      const double tau =
        is_steady(method) ?
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2)) :
          1. / std::sqrt(std::pow(sdt, 2) + std::pow(2. * u_mag / h, 2) +
                         9 * std::pow(4 * diffusivity / (h * h), 2));

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] +=
        velocity * tracer_gradient - diffusivity * tracer_laplacian;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];

          // rhs for : - D * laplacian T +  u * grad T - f=0
          local_rhs(i) -= (diffusivity * grad_phi_T_i * tracer_gradient +
                           phi_T_i * velocity * tracer_gradient -
                           scratch_data.source[q] * phi_T_i) *
                          JxW;

          local_rhs(i) -=
            tau * (strong_residual_vec[q] * (grad_phi_T_i * velocity)) * JxW;

          local_rhs(i) -=
            vdcdd *
            scalar_product(grad_phi_T_i, dcdd_factor * tracer_gradient) * JxW;
        }
    } // end loop on quadrature points
}

template class TracerAssemblerCore<2>;
template class TracerAssemblerCore<3>;

template <int dim>
void
TracerAssemblerDGCore<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &diffusivity_vector =
    scratch_data.tracer_diffusivity;

  // Loop and quadrature informations
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // assembling local matrix and right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         diffusivity = diffusivity_vector[q];
      const Tensor<1, dim> velocity    = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];
          const auto phi_T_i      = scratch_data.phi[q][i];


          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const Tensor<1, dim> grad_phi_T_j = scratch_data.grad_phi[q][j];
              const auto           phi_T_j      = scratch_data.phi[q][j];

              // Weak form : - D * laplacian T +  u * gradT + T div(u) - f=0
              // Note that the advection term has been weakened for it to appear
              // explicitly in the weak form as a boundary term.
              local_matrix(i, j) += (diffusivity * grad_phi_T_i * grad_phi_T_j -
                                     grad_phi_T_i * velocity * phi_T_j -
                                     phi_T_i * velocity_divergence) *
                                    JxW;
            }
        }
    } // end loop on quadrature points
}


template <int dim>
void
TracerAssemblerDGCore<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &diffusivity_vector =
    scratch_data.tracer_diffusivity;

  // Loop and quadrature informations
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Assembling local right hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Gather into local variables the relevant fields
      const double         diffusivity     = diffusivity_vector[q];
      const double         tracer_value    = scratch_data.tracer_values[q];
      const Tensor<1, dim> tracer_gradient = scratch_data.tracer_gradients[q];
      const Tensor<1, dim> velocity        = scratch_data.velocity_values[q];
      const double velocity_divergence = scratch_data.velocity_divergences[q];


      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_T_i      = scratch_data.phi[q][i];
          const auto grad_phi_T_i = scratch_data.grad_phi[q][i];

          // rhs for : - D * laplacian T +  u * grad T - T div(u) - f=0
          local_rhs(i) -= (diffusivity * grad_phi_T_i * tracer_gradient -
                           grad_phi_T_i * velocity * tracer_value -
                           phi_T_i * velocity_divergence * tracer_value -
                           scratch_data.source[q] * phi_T_i) *
                          JxW;
        }
    } // end loop on quadrature points
}

template class TracerAssemblerDGCore<2>;
template class TracerAssemblerDGCore<3>;


template <int dim>
void
TracerAssemblerBDF<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Loop and quadrature informations
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual = copy_data.strong_residual;
  auto &strong_jacobian = copy_data.strong_jacobian;
  auto &local_matrix    = copy_data.local_matrix;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
  std::vector<double> tracer(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      tracer[0] = scratch_data.tracer_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        tracer[p + 1] = scratch_data.previous_tracer_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += bdf_coefs[p] * tracer[p];
        }

      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          strong_jacobian[q][j] += bdf_coefs[0] * scratch_data.phi[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i = scratch_data.phi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const double phi_u_j = scratch_data.phi[q][j];

              local_matrix(i, j) += phi_u_j * phi_u_i * bdf_coefs[0] * JxW[q];
            }
        }
    }
}

template <int dim>
void
TracerAssemblerBDF<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Loop and quadrature informations
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
  std::vector<double> tracer(1 + number_of_previous_solutions(method));

  // Loop over the quadrature points
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      tracer[0] = scratch_data.tracer_values[q];
      for (unsigned int p = 0; p < number_of_previous_solutions(method); ++p)
        tracer[p + 1] = scratch_data.previous_tracer_values[p][q];

      for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
           ++p)
        {
          strong_residual[q] += (bdf_coefs[p] * tracer[p]);
        }

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const double phi_u_i     = scratch_data.phi[q][i];
          double       local_rhs_i = 0;
          for (unsigned int p = 0; p < number_of_previous_solutions(method) + 1;
               ++p)
            {
              local_rhs_i -= bdf_coefs[p] * (tracer[p] * phi_u_i);
            }
          local_rhs(i) += local_rhs_i * JxW[q];
        }
    }
}

template class TracerAssemblerBDF<2>;
template class TracerAssemblerBDF<3>;



template <int dim>
void
TracerAssemblerSIPG<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedDGMethodsCopyData  &copy_data)
{
  const double                  penalty_factor = scratch_data.penalty_factor;
  const FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values_tracer;
  const auto                   &q_points       = fe_iv.get_quadrature_points();
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
            copy_data_face.face_matrix(i, j) +=
              fe_iv.jump_in_shape_values(i, q) // [\phi_i]
              * fe_iv.shape_value((velocity_dot_n > 0.),
                                  j,
                                  q) // phi_j^{upwind}
              * velocity_dot_n       // (u . n)
              * JxW[q];              // dx

            // Assemble the diffusion term using Nitsche
            // symmetric interior penalty method. See Larson Chap. 14. P.362
            copy_data_face.face_matrix(i, j) +=
              (scratch_data.tracer_diffusivity_face[q] *
                 (-fe_iv.average_of_shape_gradients(j, q) * normals[q] *
                    fe_iv.jump_in_shape_values(i, q) -
                  fe_iv.average_of_shape_gradients(i, q) * normals[q] *
                    fe_iv.jump_in_shape_values(j, q)) +
               (scratch_data.tracer_diffusivity_face[q]) * penalty_factor *
                 fe_iv.jump_in_shape_values(j, q) *
                 fe_iv.jump_in_shape_values(i, q)) *
              JxW[q];
          }
    }
}

template <int dim>
void
TracerAssemblerSIPG<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedDGMethodsCopyData  &copy_data)
{
  const double                  penalty_factor = scratch_data.penalty_factor;
  const FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values_tracer;
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

          // Assemble the diffusion term using Nitsche symmetric interior
          // penalty method. See Larson Chap. 14. P.362
          copy_data_face.face_rhs(i) -=
            (scratch_data.tracer_diffusivity_face[q] *
               (-scratch_data.tracer_average_gradient[q] * normals[q] *
                  fe_iv.jump_in_shape_values(i, q) -
                scratch_data.tracer_value_jump[q] * normals[q] *
                  fe_iv.average_of_shape_gradients(i, q)) +
             (scratch_data.tracer_diffusivity_face[q]) * penalty_factor *
               scratch_data.tracer_value_jump[q] *
               fe_iv.jump_in_shape_values(i, q)) *
            JxW[q];
        }
    }
}

template class TracerAssemblerSIPG<2>;
template class TracerAssemblerSIPG<3>;


template <int dim>
void
TracerAssemblerBoundaryNitsche<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedDGMethodsCopyData  &copy_data)
{
  const unsigned int           boundary_index = scratch_data.boundary_index;
  const FEFaceValuesBase<dim> &fe_face =
    scratch_data.fe_interface_values_tracer.get_fe_face_values(0);
  const auto &q_points = fe_face.get_quadrature_points();

  const double beta = scratch_data.penalty_factor;

  const unsigned int         n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
  const std::vector<double> &JxW          = fe_face.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

  for (unsigned int point = 0; point < q_points.size(); ++point)
    {
      const double velocity_dot_n =
        scratch_data.face_velocity_values[point] * normals[point];
      for (unsigned int i = 0; i < n_facet_dofs; ++i)
        {
          for (unsigned int j = 0; j < n_facet_dofs; ++j)
            {
              if (boundary_conditions_tracer.type[boundary_index] ==
                  BoundaryConditions::BoundaryType::outlet)
                {
                  if (velocity_dot_n > 0)
                    copy_data.local_matrix(i, j) +=
                      fe_face.shape_value(i, point) *
                      fe_face.shape_value(j, point) * velocity_dot_n *
                      JxW[point];
                }
              else if (boundary_conditions_tracer.type[boundary_index] ==
                       BoundaryConditions::BoundaryType::tracer_dirichlet)
                {
                  if (velocity_dot_n > 0)
                    copy_data.local_matrix(i, j) +=
                      fe_face.shape_value(i, point) *
                      fe_face.shape_value(j, point) * velocity_dot_n *
                      JxW[point];

                  copy_data.local_matrix(i, j) +=
                    scratch_data.tracer_diffusivity_face[point] *
                    (-fe_face.shape_value(i, point) *
                       fe_face.shape_grad(j, point) * normals[point] -
                     fe_face.shape_value(j, point) *
                       fe_face.shape_grad(i, point) * normals[point] +
                     beta * fe_face.shape_value(i, point) *
                       fe_face.shape_value(j, point)) *
                    JxW[point];
                }
            }
        }
    }
}

template <int dim>
void
TracerAssemblerBoundaryNitsche<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedDGMethodsCopyData  &copy_data)
{
  const unsigned int           boundary_index = scratch_data.boundary_index;
  const FEFaceValuesBase<dim> &fe_face =
    scratch_data.fe_interface_values_tracer.get_fe_face_values(0);
  const auto &q_points = fe_face.get_quadrature_points();

  const double beta = scratch_data.penalty_factor;

  const unsigned int         n_facet_dofs = fe_face.get_fe().n_dofs_per_cell();
  const std::vector<double> &JxW          = fe_face.get_JxW_values();
  const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

  // If the boundary condition is an outlet, assumes that advection comes
  // out since there is no inflow
  if (boundary_conditions_tracer.type[boundary_index] ==
      BoundaryConditions::BoundaryType::outlet)
    {
      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          const double velocity_dot_n =
            scratch_data.face_velocity_values[point] * normals[point];

          for (unsigned int i = 0; i < n_facet_dofs; ++i)
            copy_data.local_rhs(i) -= fe_face.shape_value(i, point) // \phi_i
                                      * scratch_data.values_here[point] *
                                      velocity_dot_n // u . n
                                      * JxW[point];  // dx
        }
    }
  // Else the boundary condition is a dirichlet. Evaluate the function and
  // process it accordingly
  else if (boundary_conditions_tracer.type[boundary_index] ==
           BoundaryConditions::BoundaryType::tracer_dirichlet)
    {
      std::vector<double> function_value(q_points.size());
      boundary_conditions_tracer.tracer[boundary_index]->value_list(
        q_points, function_value);

      for (unsigned int point = 0; point < q_points.size(); ++point)
        {
          const double velocity_dot_n =
            scratch_data.face_velocity_values[point] * normals[point];

          for (unsigned int i = 0; i < n_facet_dofs; ++i)
            {
              if (velocity_dot_n < 0)
                copy_data.local_rhs(i) -= fe_face.shape_value(i, point) *
                                          function_value[point] *
                                          velocity_dot_n * JxW[point];

              if (velocity_dot_n > 0)
                copy_data.local_rhs(i) -= fe_face.shape_value(i, point) *
                                          scratch_data.values_here[point] *
                                          velocity_dot_n * JxW[point];

              copy_data.local_rhs(i) -=
                scratch_data.tracer_diffusivity_face[point] *
                (-fe_face.shape_value(i, point) *
                   scratch_data.gradients_here[point] * normals[point] -
                 (scratch_data.values_here[point] - function_value[point]) *
                   fe_face.shape_grad(i, point) * normals[point] +
                 beta * fe_face.shape_value(i, point) *
                   (scratch_data.values_here[point] - function_value[point])) *
                JxW[point];
            }
        }
    }
  else if (boundary_conditions_tracer.type[boundary_index] ==
             BoundaryConditions::BoundaryType::periodic ||
           boundary_conditions_tracer.type[boundary_index] ==
             BoundaryConditions::BoundaryType::periodic_neighbor)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Periodic boundary conditions cannot be used with the tracer physics when discontinous Galerkin elements are used."));
    }
  // If it is an unknown boundary condition, throw an exception. This case
  // should never occur.
  else
    {
      AssertThrow(false,
                  ExcMessage(
                    "No valid boundary conditions types were identified."));
    }
}

template class TracerAssemblerBoundaryNitsche<2>;
template class TracerAssemblerBoundaryNitsche<3>;

template <int dim>
void
TracerAssemblerReaction<dim>::assemble_matrix(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &k = scratch_data.tracer_reaction_prefactor;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  auto &strong_jacobian_vec = copy_data.strong_jacobian;
  auto &local_matrix        = copy_data.local_matrix;

  // Assembling local matrix
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double reaction_coeff = k[q];
      const double C              = scratch_data.tracer_values[q];

      // Update the strong Jacobian with the reaction term contribution:
      // It comes from the derivative -d/dC (k * C) = -(k + C*dk/dC).
      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto phi_T_j = scratch_data.phi[q][j];
          strong_jacobian_vec[q][j] +=
            -(reaction_coeff +
              scratch_data.grad_tracer_reaction_prefactor[q] * C) *
            phi_T_j;
        }

      // Store JxW in a local variable for faster access.
      const double JxW = JxW_vec[q];

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_T_i = scratch_data.phi[q][i];

          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto &phi_T_j = scratch_data.phi[q][j];

              local_matrix(i, j) += reaction_coeff * phi_T_i * phi_T_j * JxW;
            }
        }
    } // end loop on quadrature points
}

template <int dim>
void
TracerAssemblerReaction<dim>::assemble_rhs(
  const TracerScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData    &copy_data)
{
  // Scheme and physical properties
  const std::vector<double> &k = scratch_data.tracer_reaction_prefactor;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &strong_residual_vec = copy_data.strong_residual;
  auto &local_rhs           = copy_data.local_rhs;

  // Assembling right-hand side
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Store JxW in local variable for faster access
      const double JxW = JxW_vec[q];

      // Get the tracer concentration value at this quadrature point
      const double C = scratch_data.tracer_values[q];

      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] += -k[q] * C;

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto &phi_T_i = scratch_data.phi[q][i];

          // Add reaction term to the RHS
          local_rhs(i) -= k[q] * phi_T_i * C * JxW;
        }
    } // end loop on quadrature points
}

template class TracerAssemblerReaction<2>;
template class TracerAssemblerReaction<3>;
