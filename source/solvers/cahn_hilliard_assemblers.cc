#include <core/bdf.h>
#include <core/time_integration_utilities.h>


#include <solvers/copy_data.h>
#include <solvers/cahn_hilliard_assemblers.h>


template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                          StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  const double epsilon = 1.0;
  const double W = 1.0;
  const double D = 1.0;

  auto &local_matrix = copy_data.local_matrix;

  Tensor<1,dim> velocity_field;

  for (unsigned int q = 0; q<n_q_points; ++q)
    {
      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        const double phi_phase_i = scratch_data.phi_phase[q][i];
        const Tensor<1,dim> grad_phi_phase_i = scratch_data.grad_phi_phase[q][i];

        const double phi_potential_i = scratch_data.phi_potential[q][i];
        const Tensor<1,dim> grad_phi_potential_i = scratch_data.grad_phi_potential[q][i];

        const double phase_order_value = scratch_data.phase_order_values[q];

        for (unsigned int j = 0; j < n_dofs ; ++j)
          {
            const double phi_phase_j = scratch_data.phi_phase[q][j];
            const Tensor<1,dim> grad_phi_phase_j = scratch_data.grad_phi_phase[q][j];

            const double phi_potential_j = scratch_data.phi_potential[q][j];
            const Tensor<1,dim> grad_phi_potential_j = scratch_data.grad_phi_potential[q][j];

            local_matrix(i,j) +=  (velocity_field*grad_phi_phase_j*phi_phase_i //First equation
                                  + D*grad_phi_phase_i*grad_phi_potential_j
                                  + phi_potential_i*phi_potential_j
                                  - 4*W*phi_potential_i*(3*phase_order_value*phase_order_value*phi_phase_j - phi_phase_j)
                                  - epsilon*epsilon*grad_phi_potential_i*grad_phi_potential_j)*JxW;

//            std::cout<<"local matrix(i,j) = "<<local_matrix(i,j)<<std::endl;
//            std::cout<<"phi_phase_i = "<<phi_phase_i<<std::endl;
//            std::cout<<"phi_phase_j = "<<phi_phase_j<<std::endl;
//            std::cout<<"phase_order_value = "<<phase_order_value<<std::endl;
//            std::cout<<"JxW = "<<JxW<<std::endl;
//            std::cout<<"grad_phi_potential_i = "<<grad_phi_potential_i<<std::endl;
//            std::cout<<"grad_phi_potential_i = "<<grad_phi_potential_i<<std::endl;

          }
      }
    }// end loop on quadrature points
}



template <int dim>
void
CahnHilliardAssemblerCore<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                             StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  const double epsilon = 1.0;
  const double W = 1.0;
  const double D = 1.0;

  auto &local_rhs           = copy_data.local_rhs;

  Tensor<1,dim> velocity_field;

  for (unsigned int q = 0; q<n_q_points; ++q)
    {
    // Store JxW in local variable for faster access;
    const double JxW = JxW_vec[q];
    const double phase_order_value = scratch_data.phase_order_values[q];
    const Tensor<1,dim> phase_order_gradient = scratch_data.phase_order_gradients[q];
    const double potential_value = scratch_data.chemical_potential_values[q];
    const Tensor<1,dim> potential_gradient = scratch_data.chemical_potential_gradients[q];

    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        const double phi_phase_i = scratch_data.phi_phase[q][i];
        const Tensor<1,dim> grad_phi_phase_i = scratch_data.grad_phi_phase[q][i];
        const double phi_potential_i = scratch_data.phi_potential[q][i];
        const Tensor<1,dim> grad_phi_potential_i = scratch_data.grad_phi_potential[q][i];

        local_rhs(i) -=  (velocity_field*phase_order_gradient*phi_phase_i //First equation
                         + D*grad_phi_phase_i*potential_gradient
                         + phi_potential_i*potential_value
                         - 4*W*phi_potential_i*(phase_order_value*phase_order_value*phase_order_value - phase_order_value)
                         - epsilon*epsilon*grad_phi_potential_i*phase_order_gradient)*JxW;

//        std::cout<<"local rhs(i) = "<<local_rhs(i)<<std::endl;
//        std::cout<<"phi_phase_i = "<<phi_phase_i<<std::endl;
//        std::cout<<"phase_order_value = "<<phase_order_value<<std::endl;
//        std::cout<<"phase_order_gradient = "<<phase_order_gradient<<std::endl;
//        std::cout<<"potential_value = "<<potential_value<<std::endl;
//        std::cout<<"potential_gradient = "<<potential_gradient<<std::endl;
//        std::cout<<"JxW = "<<JxW<<std::endl;
//        std::cout<<"grad_phi_phase_i= "<<grad_phi_phase_i<<std::endl;
//        std::cout<<"grad_phi_potential_i = "<<grad_phi_potential_i<<std::endl;
      }
    }// end loop on quadrature points
}


template class CahnHilliardAssemblerCore<2>;
template class CahnHilliardAssemblerCore<3>;

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_matrix(CahnHilliardScratchData<dim> &scratch_data,
                                               StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
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



    for (unsigned int i = 0; i < n_dofs; ++i)
      {
        const double phi_phase_i = scratch_data.phi_phase[q][i];
        for (unsigned int j = 0; j < n_dofs; ++j)
          {
            const double phi_phase_j = scratch_data.phi_phase[q][j];

            local_matrix(i, j) += phi_phase_j * phi_phase_i * bdf_coefs[0] * JxW[q];
          }
      }
    }
}

template <int dim>
void
CahnHilliardAssemblerBDF<dim>::assemble_rhs(CahnHilliardScratchData<dim> &   scratch_data,
                                            StabilizedMethodsCopyData &copy_data)
{
  // Loop and quadrature informations
  const auto &       JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
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
