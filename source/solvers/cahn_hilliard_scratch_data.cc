// SPDX-FileCopyrightText: Copyright (c) 2023-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/cahn_hilliard_scratch_data.h>

template <int dim>
CahnHilliardScratchData<dim>::CahnHilliardScratchData(
  const PhysicalPropertiesManager &properties_manager,
  const FESystem<dim>             &fe_cahn_hilliard,
  const Quadrature<dim>           &quadrature,
  const Mapping<dim>              &mapping,
  const FiniteElement<dim>        &fe_fd,
  const Quadrature<dim - 1>       &face_quadrature)
  : properties_manager(properties_manager)
  , fe_values_cahn_hilliard(mapping,
                            fe_cahn_hilliard,
                            quadrature,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_hessians)
  , fe_values_fd(mapping, fe_fd, quadrature, update_values)
  , fe_face_values_cahn_hilliard(mapping,
                                 fe_cahn_hilliard,
                                 face_quadrature,
                                 update_values | update_quadrature_points |
                                   update_JxW_values | update_gradients |
                                   update_normal_vectors)
{
  allocate();
}

template <int dim>
CahnHilliardScratchData<dim>::CahnHilliardScratchData(
  const CahnHilliardScratchData<dim> &sd)
  : properties_manager(sd.properties_manager)
  , fe_values_cahn_hilliard(sd.fe_values_cahn_hilliard.get_mapping(),
                            sd.fe_values_cahn_hilliard.get_fe(),
                            sd.fe_values_cahn_hilliard.get_quadrature(),
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients |
                              update_hessians)
  , fe_values_fd(sd.fe_values_fd.get_mapping(),
                 sd.fe_values_fd.get_fe(),
                 sd.fe_values_fd.get_quadrature(),
                 update_values)
  , fe_face_values_cahn_hilliard(
      sd.fe_face_values_cahn_hilliard.get_mapping(),
      sd.fe_face_values_cahn_hilliard.get_fe(),
      sd.fe_face_values_cahn_hilliard.get_quadrature(),
      update_values | update_quadrature_points | update_JxW_values |
        update_gradients | update_normal_vectors)
{
  allocate();
}

template <int dim>
void
CahnHilliardScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_cahn_hilliard.get_quadrature().size();
  this->n_dofs     = fe_values_cahn_hilliard.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Forcing term array
  this->source_phase_order        = std::vector<double>(n_q_points);
  this->source_chemical_potential = std::vector<double>(n_q_points);

  // Initialize arrays related to phase order and chemical potential
  this->phase_order.component        = 0;
  this->chemical_potential.component = 1;

  // Variables for the Cahn-Hilliard equations
  this->phase_order_values            = std::vector<double>(n_q_points);
  this->phase_order_gradients         = std::vector<Tensor<1, dim>>(n_q_points);
  this->phase_order_laplacians        = std::vector<double>(n_q_points);
  this->chemical_potential_values     = std::vector<double>(n_q_points);
  this->chemical_potential_gradients  = std::vector<Tensor<1, dim>>(n_q_points);
  this->chemical_potential_laplacians = std::vector<double>(n_q_points);

  // Phase order for BDF schemes
  this->previous_phase_order_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  this->previous_chemical_potential_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Phase-order shape functions
  this->phi_phase.reinit(n_q_points, n_dofs);
  this->grad_phi_phase.reinit(n_q_points, n_dofs);
  this->hess_phi_phase.reinit(n_q_points, n_dofs);
  this->laplacian_phi_phase.reinit(n_q_points, n_dofs);

  // Chemical potential shape functions
  this->phi_potential.reinit(n_q_points, n_dofs);
  this->grad_phi_potential.reinit(n_q_points, n_dofs);
  this->hess_phi_potential.reinit(n_q_points, n_dofs);
  this->laplacian_phi_potential.reinit(n_q_points, n_dofs);

  // Velocity values
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  this->velocity_gradient_values =
    std::vector<Tensor<2, dim>>(this->n_q_points);

  // Allocate physical properties
  this->surface_tension                 = std::vector<double>(n_q_points);
  this->mobility_cahn_hilliard          = std::vector<double>(n_q_points);
  this->mobility_cahn_hilliard_gradient = std::vector<double>(n_q_points);

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));
  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));
}

template <int dim>
template <typename VectorType>
void
CahnHilliardScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions,
  Function<dim>                                        *source_function)
{
  this->phase_order.component        = 0;
  this->chemical_potential.component = 1;

  this->fe_values_cahn_hilliard.reinit(cell);

  quadrature_points = this->fe_values_cahn_hilliard.get_quadrature_points();
  auto &fe_cahn_hilliard = this->fe_values_cahn_hilliard.get_fe();

  source_function->value_list(quadrature_points, source_phase_order, 0);
  source_function->value_list(quadrature_points, source_chemical_potential, 1);

  // Compute cell diameter
  double cell_measure = compute_cell_measure_with_JxW(
    this->fe_values_cahn_hilliard.get_JxW_values());
  this->cell_size =
    compute_cell_diameter<dim>(cell_measure, fe_cahn_hilliard.degree);

  // Gather Phi and eta (values, gradient and laplacian)
  this->fe_values_cahn_hilliard[phase_order].get_function_values(
    current_solution, this->phase_order_values);
  this->fe_values_cahn_hilliard[phase_order].get_function_gradients(
    current_solution, this->phase_order_gradients);
  this->fe_values_cahn_hilliard[phase_order].get_function_laplacians(
    current_solution, this->phase_order_laplacians);

  this->fe_values_cahn_hilliard[chemical_potential].get_function_values(
    current_solution, this->chemical_potential_values);
  this->fe_values_cahn_hilliard[chemical_potential].get_function_gradients(
    current_solution, this->chemical_potential_gradients);
  this->fe_values_cahn_hilliard[chemical_potential].get_function_laplacians(
    current_solution, this->chemical_potential_laplacians);

  // Gather previous phase order values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values_cahn_hilliard[phase_order].get_function_values(
        previous_solutions[p], previous_phase_order_values[p]);
    }

  // Gather previous chemical potential values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values_cahn_hilliard[chemical_potential].get_function_values(
        previous_solutions[p], previous_chemical_potential_values[p]);
    }

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      this->JxW[q] = this->fe_values_cahn_hilliard.JxW(q);

      for (unsigned int k = 0; k < n_dofs; ++k)
        {
          // Shape functions for the phase order
          this->phi_phase[q][k] =
            this->fe_values_cahn_hilliard[phase_order].value(k, q);
          this->grad_phi_phase[q][k] =
            this->fe_values_cahn_hilliard[phase_order].gradient(k, q);
          this->hess_phi_phase[q][k] =
            this->fe_values_cahn_hilliard[phase_order].hessian(k, q);
          this->laplacian_phi_phase[q][k] = trace(this->hess_phi_phase[q][k]);

          // Shape functions for the chemical potential
          this->phi_potential[q][k] =
            this->fe_values_cahn_hilliard[chemical_potential].value(k, q);
          this->grad_phi_potential[q][k] =
            this->fe_values_cahn_hilliard[chemical_potential].gradient(k, q);
          this->hess_phi_potential[q][k] =
            this->fe_values_cahn_hilliard[chemical_potential].hessian(k, q);
          this->laplacian_phi_potential[q][k] =
            trace(this->hess_phi_potential[q][k]);
        }
    }

  this->is_boundary_cell =
    cell->at_boundary(); // The attribute needs to be updated because the
                         // assembler for the angle of contact boundary
                         // condition needs to know if the cell is at the
                         // boundary
  if (this->is_boundary_cell)
    {
      n_faces          = cell->n_faces();
      is_boundary_face = std::vector<bool>(n_faces, false);
      n_faces_q_points = fe_face_values_cahn_hilliard.get_quadrature().size();
      boundary_face_id = std::vector<unsigned int>(n_faces);

      face_JxW.reinit(n_faces, n_faces_q_points);
      this->grad_phi_face_phase.reinit(n_faces, n_faces_q_points, n_dofs);

      this->face_phase_grad_values = std::vector<std::vector<Tensor<1, dim>>>(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      this->face_normal = std::vector<std::vector<Tensor<1, dim>>>(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      for (const auto face : cell->face_indices())
        {
          this->is_boundary_face[face] = cell->face(face)->at_boundary();
          if (this->is_boundary_face[face])
            {
              fe_face_values_cahn_hilliard.reinit(cell, face);
              boundary_face_id[face] = cell->face(face)->boundary_id();

              this->fe_face_values_cahn_hilliard[phase_order]
                .get_function_gradients(current_solution,
                                        this->face_phase_grad_values[face]);

              for (unsigned int q = 0; q < n_faces_q_points; ++q)
                {
                  face_JxW[face][q] = fe_face_values_cahn_hilliard.JxW(q);
                  this->face_normal[face][q] =
                    this->fe_face_values_cahn_hilliard.normal_vector(q);
                  for (const unsigned int k :
                       fe_face_values_cahn_hilliard.dof_indices())
                    {
                      this->grad_phi_face_phase[face][q][k] =
                        this->fe_face_values_cahn_hilliard[phase_order]
                          .gradient(k, q);
                    }
                }
            }
        }
    }
}

template <int dim>
template <typename VectorType>
void
CahnHilliardScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const Parameters::ALE<dim>                           &ale)
{
  FEValuesExtractors::Vector velocities(0);
  this->fe_values_fd.reinit(cell);

  this->fe_values_fd[velocities].get_function_values(current_solution,
                                                     velocity_values);

  if (!ale.enabled())
    return;


  // ALE enabled, so extract the ALE velocity and subtract it from the
  // velocity obtained from the fluid dynamics
  Tensor<1, dim>                                  velocity_ale;
  std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
    ale.velocity;
  Vector<double> velocity_ale_vector(dim);

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      velocity_ale_function->vector_value(quadrature_points[q],
                                          velocity_ale_vector);
      for (int d = 0; d < dim; ++d)
        velocity_ale[d] = velocity_ale_vector[d];

      velocity_values[q] -= velocity_ale;
    }
}

template <int dim>
void
CahnHilliardScratchData<dim>::calculate_physical_properties()
{
  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          throw std::runtime_error("Unsupported number of fluids (<2)");
        }
      case 2:
        {
          set_field_vector(field::phase_order_cahn_hilliard,
                           this->phase_order_values,
                           this->fields);

          // Gather properties from material interactions
          const auto material_interaction_id =
            properties_manager.get_material_interaction_id(
              material_interactions_type::fluid_fluid, 0, 1);

          // Gather surface tension
          const auto surface_tension_model =
            properties_manager.get_surface_tension(material_interaction_id);
          surface_tension_model->vector_value(fields, surface_tension);

          // Gather mobility
          const auto mobility_cahn_hilliard_model =
            properties_manager.get_mobility_cahn_hilliard(
              material_interaction_id);
          mobility_cahn_hilliard_model->vector_value(fields,
                                                     mobility_cahn_hilliard);
          mobility_cahn_hilliard_model->vector_jacobian(
            fields,
            field::phase_order_cahn_hilliard,
            mobility_cahn_hilliard_gradient);

          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}

template class CahnHilliardScratchData<2>;
template class CahnHilliardScratchData<3>;

template void
CahnHilliardScratchData<2>::reinit<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions,
  Function<2>                               *source_function);

template void
CahnHilliardScratchData<3>::reinit<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions,
  Function<3>                               *source_function);

template void
CahnHilliardScratchData<2>::reinit_velocity<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const Parameters::ALE<2>                  &ale);

template void
CahnHilliardScratchData<2>::reinit_velocity<GlobalBlockVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalBlockVectorType               &current_solution,
  const Parameters::ALE<2>                  &ale);

template void
CahnHilliardScratchData<3>::reinit_velocity<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const Parameters::ALE<3>                  &ale);

template void
CahnHilliardScratchData<3>::reinit_velocity<GlobalBlockVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalBlockVectorType               &current_solution,
  const Parameters::ALE<3>                  &ale);
