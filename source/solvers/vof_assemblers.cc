#include <core/bdf.h>
#include <core/time_integration_utilities.h>

<<<<<<< HEAD
#include <solvers/free_surface.h>

#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/numerics/vector_tools.h>
=======
#include <solvers/copy_data.h>
#include <solvers/vof_assemblers.h>
>>>>>>> edc38f4 (Revise the new assembler system to fix the breaking tests)


template <int dim>
void
VOFAssemblerCore<dim>::assemble_matrix(VOFScratchData<dim> &scratch_data,
                                          StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto   method      = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
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
      const Tensor<1, dim> velocity        = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      const auto phase_gradient = scratch_data.phase_gradients[q];
      const auto phase_gradient_norm = phase_gradient.norm();

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter and the compression term for the phase
      // indicator
      const double u_mag = std::max(velocity.norm(), 1e-12);

      // Implementation of a DCDD shock capturing scheme.
      // For more information see
      // Tezduyar, T. E., & Park, Y. J. (1986). Discontinuity-capturing
      // finite element formulations for nonlinear
      // convection-diffusion-reaction equations. Computer methods in
      // applied mechanics and engineering, 59(3), 307-325.

      // Gather the order of the free surface interpolation
      const double order = this->fem_parameters.free_surface_order;

      // Calculate the artificial viscosity of the shock capture
      const double vdcdd = (0.5 * h) *
                           (velocity.norm() * velocity.norm()) *
                           pow(phase_gradient_norm * h, order);

      const double tol = 1e-12;

      // We neglect to remove the diffusion aligned with the velocity
      // as is done in the original article. We re-enable those
      // terms if artificial diffusion becomes a problem
      // Tensor<1, dim> s = velocity / (velocity.norm() + 1e-12);
      // const Tensor<2, dim> k_corr      = (r * s) * outer_product(s,
      // s);

      // Calculate the unit vector associated with the phase gradient
      Tensor<1, dim> r =
        phase_gradient / (phase_gradient_norm + tol);

      // Calculate the dyadic product of this vector with itself
      const Tensor<2, dim> rr = outer_product(r, r);
      // Agglomerate this as a factor in case we want to remove
      // the contribution aligned with the velocity
      const Tensor<2, dim> dcdd_factor = rr; // - k_corr;

      // Gradient of the shock capturing viscosity for the assemblyu
      // of the jacobian matrix
      const double d_vdcdd =
        order * (0.5 * h * h) * (velocity.norm() * velocity.norm()) *
        pow(phase_gradient_norm * h, order - 1);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step. Hypothesis : advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          h / (2. * u_mag) :
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) + std::pow(sdt, 2));


      for (unsigned int j = 0; j < n_dofs; ++j)
        {
          const auto grad_phi_phase_j = scratch_data.grad_phi[q][j];
          const auto laplacian_phi_phase_j =
            scratch_data.laplacian_phi[q][j];

          // Strong Jacobian associated with the GLS
          // stabilization
          strong_jacobian_vec[q][j] += velocity * grad_phi_phase_j;
//std::cout << strong_jacobian_vec[q][j] << std::endl;

          if (DCDD)
            strong_jacobian_vec[q][j] += -vdcdd * laplacian_phi_phase_j;

         // std::cout << strong_jacobian_vec[q][j] << std::endl;
        }



      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];


              for (unsigned int j = 0; j < n_dofs; ++j)
                {
                  const auto grad_phi_phase_j = scratch_data.grad_phi[q][j];

                  // Weak form for advection: u * nabla(phase) = 0
                  local_matrix(i, j) +=
                    (phi_phase_i * velocity * grad_phi_phase_j) * JxW;

//std::cout << (phi_phase_i * velocity * grad_phi_phase_j) * JxW << std::endl;


                 // std::cout << "paramsssssssssssssssssssssssssss    "<< tau  << "   " << strong_jacobian_vec[q][j] <<  "    " << grad_phi_phase_i <<  "    " << velocity << "   " << JxW <<std::endl;




                  // Addition to the cell matrix for GLS stabilization
                  local_matrix(i, j) +=
                    tau * strong_jacobian_vec[q][j] *
                    (grad_phi_phase_i * velocity) * JxW;

                  if (DCDD)
                    {
                      local_matrix(i, j) +=
                        (vdcdd * scalar_product(grad_phi_phase_j,
                                                dcdd_factor *
                                                  grad_phi_phase_i) +
                         d_vdcdd * grad_phi_phase_j.norm() *
                           scalar_product(phase_gradient,
                                          dcdd_factor *
                                            grad_phi_phase_i)) *
                        JxW;
                    }
                }

        }
    } // end loop on quadrature points
}



template <int dim>
void
VOFAssemblerCore<dim>::assemble_rhs(VOFScratchData<dim> &   scratch_data,
                                       StabilizedMethodsCopyData &copy_data)
{
  // Scheme and physical properties
  const auto   method      = this->simulation_control->get_assembly_method();

  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
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
      const Tensor<1, dim> phase_gradient  = scratch_data.phase_gradients[q];
      const double phase_laplacians = scratch_data.phase_laplacians[q];
      const Tensor<1, dim> velocity         = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      // Shock capturing viscosity term
      const double order = scratch_data.fe_values_fs.get_fe().degree;


      const double vdcdd = (0.5 * h) * (velocity.norm() * velocity.norm()) *
                           pow(phase_gradient.norm() * h, order);

      const double   tolerance = 1e-12;
      Tensor<1, dim> s         = velocity / (velocity.norm() + tolerance);
      Tensor<1, dim> r = phase_gradient / (phase_gradient.norm() + tolerance);

      const Tensor<2, dim> k_corr      = (r * s) * outer_product(s, s);
      const Tensor<2, dim> rr          = outer_product(r, r);
      const Tensor<2, dim> dcdd_factor = rr - k_corr;

      // Calculation of the magnitude of the velocity for the
      // stabilization parameter
      const double u_mag = std::max(velocity.norm(), tolerance);

      // Calculation of the GLS stabilization parameter. The
      // stabilization parameter used is different if the simulation is
      // steady or unsteady. In the unsteady case it includes the value
      // of the time-step. Hypothesis : advection dominated problem
      // (Pe>3) [Bochev et al., Stability of the SUPG finite element
      // method for transient advection-diffusion problems, CMAME 2004]
      const double tau =
        is_steady(method) ?
          h / (2. * u_mag) :
          1. / std::sqrt(std::pow(2. * u_mag / h, 2) + std::pow(sdt, 2));







      // Calculate the strong residual for GLS stabilization
      strong_residual_vec[q] +=
        velocity * phase_gradient;

      if (DCDD)
        strong_residual_vec[q] += -vdcdd * phase_laplacians;


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          const auto grad_phi_phase_i = scratch_data.grad_phi[q][i];


          // rhs for : u * nabla(phase) = 0
          local_rhs(i) -=
            (phi_phase_i * velocity * phase_gradient) *
            JxW;


          local_rhs(i) -= tau *
                         (strong_residual_vec[q] *
                          (grad_phi_phase_i * velocity)) *
                         JxW;

          if (DCDD)
            {
              local_rhs(i) +=
                -vdcdd *
                scalar_product(phase_gradient,
                               dcdd_factor * grad_phi_phase_i) *
                JxW;
            }

      }

  }
}

<<<<<<< HEAD
template <int dim>
void
FreeSurface<dim>::compute_kelly(dealii::Vector<float> &estimated_error_per_cell)
{
  if (this->simulation_parameters.mesh_adaptation.variable ==
      Parameters::MeshAdaptation::Variable::phase)
    {
      const FEValuesExtractors::Scalar phase(0);

      KellyErrorEstimator<dim>::estimate(
        *this->mapping,
        dof_handler,
        *this->face_quadrature,
        typename std::map<types::boundary_id, const Function<dim, double> *>(),
        present_solution,
        estimated_error_per_cell,
        fe->component_mask(phase));
    }
}

template <int dim>
void
FreeSurface<dim>::write_checkpoint()
{
  std::vector<const TrilinosWrappers::MPI::Vector *> sol_set_transfer;
=======
template class VOFAssemblerCore<2>;
template class VOFAssemblerCore<3>;

>>>>>>> edc38f4 (Revise the new assembler system to fix the breaking tests)


template <int dim>
void
VOFAssemblerBDF<dim>::assemble_matrix(VOFScratchData<dim> &scratch_data,
                                         StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
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

          strong_jacobian[q][j] +=
                   bdf_coefs[0] * scratch_data.phi[q][j];
        }


      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          const auto phi_phase_i      = scratch_data.phi[q][i];
          for (unsigned int j = 0; j < n_dofs; ++j)
            {
              const auto phi_phase_j      = scratch_data.phi[q][j];

              local_matrix(i, j) += phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];



            }
        }
  }
<<<<<<< HEAD
  zero_constraints.close();

  // Sparse matrices initialization
  DynamicSparsityPattern dsp(this->dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(this->dof_handler,
                                  dsp,
                                  nonzero_constraints,
                                  /*keep_constrained_dofs = */ true);

  SparsityTools::distribute_sparsity_pattern(dsp,
                                             locally_owned_dofs,
                                             mpi_communicator,
                                             locally_relevant_dofs);
  system_matrix.reinit(locally_owned_dofs,
                       locally_owned_dofs,
                       dsp,
                       mpi_communicator);

  this->pcout << "   Number of free surface degrees of freedom: "
              << dof_handler.n_dofs() << std::endl;

  // Provide the free surface dof_handler and solution pointers to the
  // multiphysics interface
  multiphysics->set_dof_handler(PhysicsID::free_surface, &dof_handler);
  multiphysics->set_solution(PhysicsID::free_surface, &present_solution);
  // the fluid at present iteration is solved before the free surface, and
  // after percolate is called for the previous iteration, so at the time the
  // getter is called solution_m2 = solution_m1
  // TODO deactivated for now (inertia is considered with a constant density),
  // see if needed / to be debugged
  multiphysics->set_solution_m1(PhysicsID::free_surface,
                                &previous_solutions[0]);
=======
>>>>>>> edc38f4 (Revise the new assembler system to fix the breaking tests)
}




template <int dim>
void
VOFAssemblerBDF<dim>::assemble_rhs(VOFScratchData<dim> &   scratch_data,
                                      StabilizedMethodsCopyData & copy_data)
{
    // Loop and quadrature informations
    const auto &       JxW        = scratch_data.JxW;
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
            const double phi_phase_i     = scratch_data.phi[q][i];
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
