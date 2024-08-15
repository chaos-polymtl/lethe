/*---------------------------------------------------------------------
 *
 * Copyright (C) 2019 - by the Lethe authors
 *
 * This file is part of the Lethe library
 *
 * The Lethe library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the Lethe distribution.
 *
 * ---------------------------------------------------------------------
 */

#include <core/bdf.h>
#include <core/time_integration_utilities.h>

#include <solvers/copy_data.h>
#include <solvers/reactive_species_assemblers.h>

template <int dim>
void
ReactiveSpeciesAssemblerCore<dim>::assemble_matrix(
  ReactiveSpeciesScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData       &copy_data)
{
  // Gather physical properties
  const double epsilon   = this->epsilon;
  const double cell_size = scratch_data.cell_size;
  const double xi =
    this->reactive_species_parameters.potential_smoothing_coefficient;

  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double lambda =
        3 * epsilon * scratch_data.surface_tension[q] / (2 * sqrt(2));

      const double mobility = scratch_data.mobility_reactive_species[q];
      const double mobility_gradient =
        scratch_data.mobility_reactive_species_gradient[q];

      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW = JxW_vec[q];

      const double phase_order_value = scratch_data.phase_order_values[q];
      const Tensor<1, dim> potential_gradient =
        scratch_data.chemical_potential_gradients[q];

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
              const double         phi_phase_j = scratch_data.phi_phase[q][j];
              const Tensor<1, dim> grad_phi_phase_j =
                scratch_data.grad_phi_phase[q][j];

              const double phi_potential_j = scratch_data.phi_potential[q][j];
              const Tensor<1, dim> grad_phi_potential_j =
                scratch_data.grad_phi_potential[q][j];

              local_matrix(i, j) +=
                // Phase order equation
                (phi_phase_i * (velocity_field * grad_phi_phase_j) +
                 mobility * grad_phi_phase_i * grad_phi_potential_j +
                 mobility_gradient * grad_phi_phase_i * potential_gradient *
                   phi_phase_j
                 // Chemical potential equation
                 + phi_potential_i * phi_potential_j -
                 (lambda / (epsilon * epsilon)) * phi_potential_i *
                   (3 * phase_order_value * phase_order_value - 1.0) *
                   phi_phase_j -
                 lambda * grad_phi_potential_i * grad_phi_phase_j
                 // Chemical potential smoothing
                 + xi * cell_size * cell_size * grad_phi_potential_i *
                     grad_phi_potential_j) *
                JxW;
            }
        }
    }
}



template <int dim>
void
ReactiveSpeciesAssemblerCore<dim>::assemble_rhs(
  ReactiveSpeciesScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData       &copy_data)
{
  // Gather physical properties
  const double epsilon   = this->epsilon;
  const double cell_size = scratch_data.cell_size;
  const double xi =
    this->reactive_species_parameters.potential_smoothing_coefficient;


  // Loop and quadrature information
  const auto        &JxW_vec    = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      const double lambda =
        3 * epsilon * scratch_data.surface_tension[q] / (2 * sqrt(2));
      const double mobility = scratch_data.mobility_reactive_species[q];


      // Gather into local variables the relevant fields
      const Tensor<1, dim> velocity_field = scratch_data.velocity_values[q];

      // Store JxW in local variable for faster access;
      const double JxW               = JxW_vec[q];
      const double phase_order_value = scratch_data.phase_order_values[q];
      const Tensor<1, dim> phase_order_gradient =
        scratch_data.phase_order_gradients[q];
      const double potential_value = scratch_data.chemical_potential_values[q];
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
            // Phase order equation
            (-phi_phase_i * (velocity_field * phase_order_gradient) -
             mobility * grad_phi_phase_i * potential_gradient
             // Chemical potential equation
             - phi_potential_i * potential_value +
             (lambda / (epsilon * epsilon)) * phi_potential_i *
               (phase_order_value * phase_order_value - 1) * phase_order_value +
             lambda * grad_phi_potential_i * phase_order_gradient
             // Chemical potential smoothing
             - xi * cell_size * cell_size * grad_phi_potential_i *
                 potential_gradient
             // Source term
             + source_phase_order * phi_phase_i +
             source_chemical_potential * phi_potential_i) *
            JxW;
        }
    }
}

template class ReactiveSpeciesAssemblerCore<2>;
template class ReactiveSpeciesAssemblerCore<3>;

template <int dim>
void
ReactiveSpeciesAssemblerBDF<dim>::assemble_matrix(
  ReactiveSpeciesScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData       &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_matrix = copy_data.local_matrix;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
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
ReactiveSpeciesAssemblerBDF<dim>::assemble_rhs(
  ReactiveSpeciesScratchData<dim> &scratch_data,
  StabilizedMethodsCopyData       &copy_data)
{
  // Loop and quadrature information
  const auto        &JxW        = scratch_data.JxW;
  const unsigned int n_q_points = scratch_data.n_q_points;
  const unsigned int n_dofs     = scratch_data.n_dofs;

  // Copy data elements
  auto &local_rhs = copy_data.local_rhs;

  // Time stepping information
  const auto method = this->simulation_control->get_assembly_method();
  // Vector for the BDF coefficients
  const Vector<double> &bdf_coefs =
    this->simulation_control->get_bdf_coefficients();
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

template class ReactiveSpeciesAssemblerBDF<2>;
template class ReactiveSpeciesAssemblerBDF<3>;
