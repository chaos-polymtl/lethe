// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/navier_stokes_scratch_data.h>

#include <cfloat>

template <int dim>
NavierStokesScratchData<dim>::NavierStokesScratchData(
  const std::shared_ptr<SimulationControl> &simulation_control,
  const PhysicalPropertiesManager          &properties_manager,
  const FESystem<dim>                      &fe,
  const Quadrature<dim>                    &quadrature,
  const Mapping<dim>                       &mapping,
  const Quadrature<dim - 1>                &face_quadrature)
  : simulation_control(simulation_control)
  , properties_manager(properties_manager)
  , fe_values(mapping,
              fe,
              quadrature,
              update_values | update_quadrature_points | update_JxW_values |
                update_gradients | update_hessians | update_inverse_jacobians)
  , fe_face_values(mapping,
                   fe,
                   face_quadrature,
                   update_values | update_quadrature_points |
                     update_JxW_values | update_gradients | update_hessians |
                     update_normal_vectors)
  , sdirk_table(
      time_stepping_is_sdirk(simulation_control->get_assembly_method()) ?
        ::sdirk_table(simulation_control->get_assembly_method()) :
        SDIRKTable())
{
  allocate();

  // By default, the assembly of variables belonging to auxiliary physics is
  // disabled.
  gather_cls                                = false;
  gather_projected_phase_indicator_gradient = false;
  gather_curvature                          = false;
  gather_void_fraction                      = false;
  gather_particles_information              = false;
  gather_temperature                        = false;
  gather_cahn_hilliard                      = false;
  gather_ale                                = false;
  gather_mortar                             = false;
  gather_particle_field_project             = false;
  gather_hessian = properties_manager.is_non_newtonian();
}

template <int dim>
NavierStokesScratchData<dim>::NavierStokesScratchData(
  const NavierStokesScratchData<dim> &sd)
  : simulation_control(sd.simulation_control)
  , properties_manager(sd.properties_manager)
  , fe_values(sd.fe_values.get_mapping(),
              sd.fe_values.get_fe(),
              sd.fe_values.get_quadrature(),
              update_values | update_quadrature_points | update_JxW_values |
                update_gradients | update_hessians | update_inverse_jacobians)
  , fe_face_values(sd.fe_face_values.get_mapping(),
                   sd.fe_face_values.get_fe(),
                   sd.fe_face_values.get_quadrature(),
                   update_values | update_quadrature_points |
                     update_JxW_values | update_gradients | update_hessians |
                     update_normal_vectors)
  , sdirk_table(sd.sdirk_table)
{
  allocate();

  // By default, the assembly of variables belonging to auxiliary physics is
  // disabled.
  gather_cls                                = false;
  gather_projected_phase_indicator_gradient = false;
  gather_curvature                          = false;
  gather_void_fraction                      = false;
  gather_particles_information              = false;
  gather_temperature                        = false;
  gather_cahn_hilliard                      = false;
  gather_ale                                = false;
  gather_mortar                             = false;
  gather_particle_field_project             = false;
  gather_hessian = properties_manager.is_non_newtonian();

  if (sd.gather_cls)
    enable_cls(sd.fe_values_cls->get_fe(),
               sd.fe_values_cls->get_quadrature(),
               sd.fe_values_cls->get_mapping(),
               sd.filter);
  if (sd.gather_projected_phase_indicator_gradient)
    enable_projected_phase_indicator_gradient(
      sd.fe_values_projected_phase_indicator_gradient->get_fe(),
      sd.fe_values_projected_phase_indicator_gradient->get_quadrature(),
      sd.fe_values_projected_phase_indicator_gradient->get_mapping());
  if (sd.gather_curvature)
    enable_curvature(sd.fe_values_curvature->get_fe(),
                     sd.fe_values_curvature->get_quadrature(),
                     sd.fe_values_curvature->get_mapping());

  if (sd.gather_void_fraction)
    enable_void_fraction(sd.fe_values_void_fraction->get_fe(),
                         sd.fe_values_void_fraction->get_quadrature(),
                         sd.fe_values_void_fraction->get_mapping());
  if (sd.gather_particles_information)
    enable_particle_fluid_interactions(sd.max_number_of_particles_per_cell,
                                       sd.interpolated_void_fraction);
  if (sd.gather_temperature)
    enable_heat_transfer(sd.fe_values_temperature->get_fe(),
                         sd.fe_values_temperature->get_quadrature(),
                         sd.fe_values_temperature->get_mapping());
  if (sd.gather_cahn_hilliard)
    enable_cahn_hilliard(sd.fe_values_cahn_hilliard->get_fe(),
                         sd.fe_values_cahn_hilliard->get_quadrature(),
                         sd.fe_values_cahn_hilliard->get_mapping(),
                         sd.cahn_hilliard_filter);
  if (sd.gather_ale)
    enable_ale();
  if (sd.gather_mortar)
    enable_mortar();

  if (sd.gather_particle_field_project)
    enable_particle_field_projection(
      sd.fe_values.get_quadrature(),
      sd.fe_values.get_mapping(),
      sd.fe_values_particle_drag->get_fe(),
      sd.fe_values_particle_two_way_coupling_force->get_fe(),
      sd.fe_values_particle_velocity->get_fe(),
      sd.fe_values_particle_momentum_transfer_coefficient->get_fe());

  gather_hessian = sd.gather_hessian;
}

template <int dim>
void
NavierStokesScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values.get_quadrature().size();
  this->n_dofs     = fe_values.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Initialize component array
  this->components = std::vector<unsigned int>(n_dofs);

  // Forcing term array
  this->rhs_force =
    std::vector<Vector<double>>(n_q_points, Vector<double>(dim + 1));
  this->force       = std::vector<Tensor<1, dim>>(n_q_points);
  this->mass_source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities.first_vector_component = 0;
  this->pressure.component                = dim;

  // Velocity
  this->velocity_values              = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences         = std::vector<double>(n_q_points);
  this->velocity_gradients           = std::vector<Tensor<2, dim>>(n_q_points);
  this->velocity_laplacians          = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_hessians            = std::vector<Tensor<3, dim>>(n_q_points);
  this->velocity_gradient_divergence = std::vector<Tensor<1, dim>>(n_q_points);
  this->shear_rate                   = std::vector<double>(n_q_points);
  this->advective_velocity           = std::vector<Tensor<1, dim>>(n_q_points);

  // For SDIRK method: sum(a_ij * k_j)
  if (this->simulation_control->is_sdirk())
    {
      this->sdirk_stage_sum = std::vector<Tensor<1, dim>>(n_q_points);
    }

  // Velocity for BDF schemes
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(n_q_points));

  // Pressure
  this->pressure_values         = std::vector<double>(n_q_points);
  this->pressure_gradients      = std::vector<Tensor<1, dim>>(n_q_points);
  this->pressure_scaling_factor = 1;

  // Pressure for BDF schemes (compressible Navier-Stokes)
  this->previous_pressure_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));
  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi_u.reinit(n_q_points, n_dofs);
  this->grad_phi_u.reinit(n_q_points, n_dofs);
  this->div_phi_u.reinit(n_q_points, n_dofs);
  this->hess_phi_u.reinit(n_q_points, n_dofs);
  this->laplacian_phi_u.reinit(n_q_points, n_dofs);
  this->gradient_divergence_phi_u.reinit(n_q_points, n_dofs);

  // Pressure shape functions
  this->phi_p.reinit(n_q_points, n_dofs);
  this->grad_phi_p.reinit(n_q_points, n_dofs);

  // Physical properties
  fields.insert(
    std::pair<field, std::vector<double>>(field::pressure, n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::shear_rate, n_q_points));

  density                               = std::vector<double>(n_q_points);
  dynamic_viscosity                     = std::vector<double>(n_q_points);
  kinematic_viscosity                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization   = std::vector<double>(n_q_points);
  kinematic_viscosity_for_stabilization = std::vector<double>(n_q_points);
  thermal_expansion                     = std::vector<double>(n_q_points);
  density_thermal_expansion             = std::vector<double>(n_q_points);
  grad_kinematic_viscosity_shear_rate   = std::vector<double>(n_q_points);

  previous_density =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
}


template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions,
  const VectorType                     &sum_over_previous_stages,
  const std::shared_ptr<Function<dim>> &forcing_function,
  const Tensor<1, dim>                 &beta_force,
  const double                          pressure_scaling_factor)
{
  this->fe_values.reinit(cell);

  quadrature_points = this->fe_values.get_quadrature_points();
  auto &fe          = this->fe_values.get_fe();

  forcing_function->vector_value_list(quadrature_points, this->rhs_force);

  // Establish the force vector
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (int d = 0; d < dim; ++d)
        {
          const unsigned int component_i =
            fe.system_to_component_index(d).first;
          this->force[q][d] = this->rhs_force[q](component_i);
        }

      const unsigned int component_mass =
        fe.system_to_component_index(dim).first;
      this->mass_source[q] = this->rhs_force[q](component_mass);

      // Correct force to include the dynamic forcing term for flow
      // control
      force[q] = force[q] + beta_force;
    }

  for (const unsigned int k : fe_values.dof_indices())
    {
      components[k] = fe.system_to_component_index(k).first;
    }

  // Compute cell diameter
  double cell_measure =
    compute_cell_measure_with_JxW(this->fe_values.get_JxW_values());
  this->cell_size = compute_cell_diameter<dim>(cell_measure, fe.degree);

  // For the SDIRK methods, \sum_{j=1}^{i-1} a_{ij} k_j is needed for the
  // assembler
  if (this->simulation_control->is_sdirk())
    {
      this->fe_values[velocities].get_function_values(sum_over_previous_stages,
                                                      this->sdirk_stage_sum);
    }

  // Gather velocity (values, gradient and laplacian)
  this->fe_values[velocities].get_function_values(current_solution,
                                                  this->velocity_values);
  this->fe_values[velocities].get_function_gradients(current_solution,
                                                     this->velocity_gradients);
  this->fe_values[velocities].get_function_laplacians(
    current_solution, this->velocity_laplacians);
  if (gather_hessian)
    this->fe_values[velocities].get_function_hessians(current_solution,
                                                      this->velocity_hessians);

  // Gather convective velocity (same as velocity_values unless ALE is
  // enabled. The ALE correction is made within the correponding reinit
  // function)
  this->advective_velocity = this->velocity_values;

  for (unsigned int q = 0; q < this->n_q_points; ++q)
    {
      this->velocity_divergences[q] = trace(this->velocity_gradients[q]);
      // In the VANS momentum equations, the viscous stress tensor includes
      // the gradient of the divergence of the velocity, ∇(∇·u) = ∂i(∂juj). We
      // therefore compute this quantity here when assembling VANS-specific
      // contributions (gated by gather_void_fraction).
      if (gather_void_fraction)
        {
          // Compute term ∂j(∂iuj) = ∂i(∂juj) for use in the strong residual
          // of the VANS momentum equations.
          this->velocity_gradient_divergence[q] = 0;
          for (int d1 = 0; d1 < dim; ++d1)
            {
              for (int d2 = 0; d2 < dim; ++d2)
                {
                  this->velocity_gradient_divergence[q][d1] +=
                    (this->velocity_hessians[q][d2][d1][d2]);
                }
            }
        }
    }

  // Gather pressure (values, gradient)
  fe_values[pressure].get_function_values(current_solution,
                                          this->pressure_values);
  fe_values[pressure].get_function_gradients(current_solution,
                                             this->pressure_gradients);
  this->pressure_scaling_factor = pressure_scaling_factor;

  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values[velocities].get_function_values(
        previous_solutions[p], previous_velocity_values[p]);
    }

  // Only gather the pressure when a pressure history is necessary
  // (compressible Navier-Stokes)
  if (!this->properties_manager.density_is_constant())
    for (unsigned int p = 0; p < previous_solutions.size(); ++p)
      {
        this->fe_values[pressure].get_function_values(
          previous_solutions[p], previous_pressure_values[p]);
      }

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      this->JxW[q] = this->fe_values.JxW(q);
      for (unsigned int k = 0; k < n_dofs; ++k)
        {
          // Velocity
          this->phi_u[q][k]      = this->fe_values[velocities].value(k, q);
          this->div_phi_u[q][k]  = this->fe_values[velocities].divergence(k, q);
          this->grad_phi_u[q][k] = this->fe_values[velocities].gradient(k, q);
          this->hess_phi_u[q][k] = this->fe_values[velocities].hessian(k, q);
          for (int d = 0; d < dim; ++d)
            this->laplacian_phi_u[q][k][d] = trace(this->hess_phi_u[q][k][d]);

          // In the VANS momentum equations, the viscous stress tensor
          // includes the gradient of the divergence of the velocity, ∇(∇·u) =
          // ∂i(∂juj). This results in a term ∂i(∂jδuj), which we compute here
          // when assembling VANS-specific contributions (gated by
          // gather_void_fraction).
          if (gather_void_fraction)
            {
              this->gradient_divergence_phi_u[q][k] = 0;
              // Compute term ∂j(∂iδuj) = ∂i(∂jδuj) to be used in the strong
              // Jacobian of the VANS momentum equations

              for (int d1 = 0; d1 < dim; ++d1)
                {
                  for (int d2 = 0; d2 < dim; ++d2)
                    {
                      this->gradient_divergence_phi_u[q][k][d1] +=
                        (this->hess_phi_u[q][k][d2][d1][d2]);
                    }
                }
            }
          // Pressure
          this->phi_p[q][k]      = this->fe_values[pressure].value(k, q);
          this->grad_phi_p[q][k] = this->fe_values[pressure].gradient(k, q);
        }
    }

  is_boundary_cell = cell->at_boundary();
  if (is_boundary_cell)
    {
      n_faces          = cell->n_faces();
      is_boundary_face = std::vector<bool>(n_faces, false);
      n_faces_q_points = fe_face_values.get_quadrature().size();
      boundary_face_id = std::vector<unsigned int>(n_faces);

      face_JxW.reinit(n_faces, n_faces_q_points);

      // Velocity and pressure values
      // First vector is face number, second quadrature point
      this->face_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      this->face_velocity_divergences =
        std::vector<std::vector<double>>(n_faces,
                                         std::vector<double>(n_faces_q_points));

      this->face_velocity_gradients = std::vector<std::vector<Tensor<2, dim>>>(
        n_faces, std::vector<Tensor<2, dim>>(n_faces_q_points));

      this->face_velocity_laplacians = std::vector<std::vector<Tensor<1, dim>>>(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      this->face_pressure_values =
        std::vector<std::vector<double>>(n_faces,
                                         std::vector<double>(n_faces_q_points));

      this->face_pressure_gradients = std::vector<std::vector<Tensor<1, dim>>>(
        n_faces, std::vector<Tensor<1, dim>>(n_faces_q_points));

      this->face_normal.reinit(n_faces, n_faces_q_points);

      this->face_quadrature_points = std::vector<std::vector<Point<dim>>>(
        n_faces, std::vector<Point<dim>>(n_faces_q_points));

      this->face_div_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

      this->face_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);
      this->face_hess_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

      this->face_laplacian_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);

      this->face_grad_phi_u.reinit(n_faces, n_faces_q_points, n_dofs);
      this->face_phi_p.reinit(n_faces, n_faces_q_points, n_dofs);
      this->face_grad_phi_p.reinit(n_faces, n_faces_q_points, n_dofs);

      for (const auto face : cell->face_indices())
        {
          is_boundary_face[face] = cell->face(face)->at_boundary();
          if (is_boundary_face[face])
            {
              fe_face_values.reinit(cell, face);
              n_dofs = fe_face_values.get_fe().n_dofs_per_cell();
              boundary_face_id[face] = cell->face(face)->boundary_id();
              // Shape functions
              // First vector is face number, second quadrature point, third
              // DOF

              // Gather velocity (values, gradient and laplacian)
              this->fe_face_values[velocities].get_function_values(
                current_solution, this->face_velocity_values[face]);
              this->fe_face_values[velocities].get_function_gradients(
                current_solution, this->face_velocity_gradients[face]);
              this->fe_face_values[velocities].get_function_laplacians(
                current_solution, this->face_velocity_laplacians[face]);

              // Gather pressure (values, gradient)
              this->fe_face_values[pressure].get_function_values(
                current_solution, this->face_pressure_values[face]);
              this->fe_face_values[pressure].get_function_gradients(
                current_solution, this->face_pressure_gradients[face]);

              for (unsigned int q = 0; q < n_faces_q_points; ++q)
                {
                  this->face_JxW[face][q] = this->fe_face_values.JxW(q);
                  this->face_normal[face][q] =
                    this->fe_face_values.normal_vector(q);
                  this->face_quadrature_points[face][q] =
                    this->fe_face_values.quadrature_point(q);
                  for (const unsigned int k : fe_face_values.dof_indices())
                    {
                      // Velocity
                      this->face_phi_u[face][q][k] =
                        this->fe_face_values[velocities].value(k, q);
                      this->face_div_phi_u[face][q][k] =
                        this->fe_face_values[velocities].divergence(k, q);
                      this->face_grad_phi_u[face][q][k] =
                        this->fe_face_values[velocities].gradient(k, q);
                      this->face_hess_phi_u[face][q][k] =
                        this->fe_face_values[velocities].hessian(k, q);
                      for (int d = 0; d < dim; ++d)
                        this->face_laplacian_phi_u[face][q][k][d] =
                          trace(this->face_hess_phi_u[face][q][k][d]);
                      // Pressure
                      this->face_phi_p[face][q][k] =
                        this->fe_face_values[pressure].value(k, q);
                      this->face_grad_phi_p[face][q][k] =
                        this->fe_face_values[pressure].gradient(k, q);
                    }
                }
            }
        }
    }
}



template <int dim>
void
NavierStokesScratchData<dim>::enable_cls(
  const FiniteElement<dim>          &fe,
  const Quadrature<dim>             &quadrature,
  const Mapping<dim>                &mapping,
  const Parameters::CLS_PhaseFilter &phase_filter_parameters)
{
  gather_cls    = true;
  fe_values_cls = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);
  point_evaluator_cls =
    std::make_shared<FEPointEvaluation<1, dim>>(mapping, fe, update_values);

  // Allocate CLS values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);
  compressibility_multiplier            = std::vector<double>(n_q_points);

  // Create filter
  filter = ConservativeLevelSetFilterBase::model_cast(phase_filter_parameters);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_cls(
  const FiniteElement<dim>                              &fe,
  const Quadrature<dim>                                 &quadrature,
  const Mapping<dim>                                    &mapping,
  const std::shared_ptr<ConservativeLevelSetFilterBase> &filter)
{
  gather_cls    = true;
  fe_values_cls = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);
  point_evaluator_cls =
    std::make_shared<FEPointEvaluation<1, dim>>(mapping, fe, update_values);

  // Allocate CLS values
  phase_values          = std::vector<double>(this->n_q_points);
  filtered_phase_values = std::vector<double>(this->n_q_points);
  previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  // For STF calculation
  phase_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  filtered_phase_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);
  compressibility_multiplier            = std::vector<double>(n_q_points);

  // Create filter
  this->filter = filter;
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_cahn_hilliard(
  const FiniteElement<dim>       &fe,
  const Quadrature<dim>          &quadrature,
  const Mapping<dim>             &mapping,
  const Parameters::CahnHilliard &cahn_hilliard_parameters)
{
  gather_cahn_hilliard    = true;
  fe_values_cahn_hilliard = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CahnHilliard values
  phase_order_cahn_hilliard_values = std::vector<double>(this->n_q_points);
  chemical_potential_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  // Allocate CahnHilliard gradients
  phase_order_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // For STF calculation
  filtered_phase_order_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));
  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);

  // Create filter
  cahn_hilliard_filter =
    CahnHilliardFilterBase::model_cast(cahn_hilliard_parameters);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_cahn_hilliard(
  const FiniteElement<dim>                      &fe,
  const Quadrature<dim>                         &quadrature,
  const Mapping<dim>                            &mapping,
  const std::shared_ptr<CahnHilliardFilterBase> &cahn_hilliard_filter)
{
  gather_cahn_hilliard    = true;
  fe_values_cahn_hilliard = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  // Allocate CahnHilliard values
  phase_order_cahn_hilliard_values = std::vector<double>(this->n_q_points);
  chemical_potential_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  // Allocate CahnHilliard gradients
  phase_order_cahn_hilliard_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  // For STF calculation
  filtered_phase_order_cahn_hilliard_values =
    std::vector<double>(this->n_q_points);

  fields.insert(
    std::pair<field, std::vector<double>>(field::phase_order_cahn_hilliard,
                                          n_q_points));

  fields.insert(std::pair<field, std::vector<double>>(
    field::phase_order_cahn_hilliard_filtered, n_q_points));

  // Allocate physical properties
  density_0                             = std::vector<double>(n_q_points);
  density_1                             = std::vector<double>(n_q_points);
  dynamic_viscosity_0                   = std::vector<double>(n_q_points);
  dynamic_viscosity_1                   = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_0 = std::vector<double>(n_q_points);
  dynamic_viscosity_for_stabilization_1 = std::vector<double>(n_q_points);
  thermal_expansion_0                   = std::vector<double>(n_q_points);
  thermal_expansion_1                   = std::vector<double>(n_q_points);
  surface_tension                       = std::vector<double>(n_q_points);
  surface_tension_gradient              = std::vector<double>(n_q_points);

  // Create filter
  this->cahn_hilliard_filter = cahn_hilliard_filter;
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_cahn_hilliard(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const VectorType &current_filtered_solution)
{
  this->fe_values_cahn_hilliard->reinit(cell);
  this->phase_order.component        = 0;
  this->chemical_potential.component = 1;

  // Gather phase indicator (values, gradients)
  this->fe_values_cahn_hilliard->operator[](phase_order)
    .get_function_values(current_solution,
                         this->phase_order_cahn_hilliard_values);
  this->fe_values_cahn_hilliard->operator[](chemical_potential)
    .get_function_values(current_solution,
                         this->chemical_potential_cahn_hilliard_values);
  this->fe_values_cahn_hilliard->operator[](phase_order)
    .get_function_gradients(current_solution,
                            this->phase_order_cahn_hilliard_gradients);

  // Gather filtered CLS solution (values, gradients)
  this->fe_values_cahn_hilliard->operator[](phase_order)
    .get_function_values(current_filtered_solution,
                         this->filtered_phase_order_cahn_hilliard_values);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_projected_phase_indicator_gradient(
  const FiniteElement<dim> &fe_projected_phase_indicator_gradient,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_projected_phase_indicator_gradient = true;
  fe_values_projected_phase_indicator_gradient =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_projected_phase_indicator_gradient,
                                    quadrature,
                                    update_values | update_gradients);

  // phase indicator gradient (PIG)
  projected_phase_indicator_gradient_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_curvature(
  const FiniteElement<dim> &fe_curvature,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_curvature    = true;
  fe_values_curvature = std::make_shared<FEValues<dim>>(
    mapping, fe_curvature, quadrature, update_values | update_gradients);

  // curvature
  curvature_values = std::vector<double>(this->n_q_points);
}


template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_cls(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const VectorType              &current_filtered_solution,
  const std::vector<VectorType> &previous_solutions)
{
  Assert(
    gather_cls,
    ExcMessage(
      "You are trying to reinit the CLS model in a cell, but you did not enable CLS for the scratch data (gather_cls=false)."));

  this->fe_values_cls->reinit(cell);
  // Gather phase indicator (values, gradient)
  this->fe_values_cls->get_function_values(current_solution,
                                           this->phase_values);
  this->fe_values_cls->get_function_values(current_filtered_solution,
                                           this->filtered_phase_values);
  this->fe_values_cls->get_function_gradients(
    current_filtered_solution, this->filtered_phase_gradient_values);
  this->fe_values_cls->get_function_gradients(current_solution,
                                              this->phase_gradient_values);

  // Gather previous phase indicator values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values_cls->get_function_values(previous_solutions[p],
                                               previous_phase_values[p]);
    }
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_projected_phase_indicator_gradient(
  const typename DoFHandler<dim>::active_cell_iterator
                   &projected_phase_indicator_gradient_cell,
  const VectorType &current_projected_phase_indicator_gradient_solution)
{
  this->fe_values_projected_phase_indicator_gradient->reinit(
    projected_phase_indicator_gradient_cell);

  FEValuesExtractors::Vector pig(0);
  // Gather phase indicator gradient
  (*fe_values_projected_phase_indicator_gradient)[pig].get_function_values(
    current_projected_phase_indicator_gradient_solution,
    this->projected_phase_indicator_gradient_values);
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_curvature(
  const typename DoFHandler<dim>::active_cell_iterator &curvature_cell,
  const VectorType &current_curvature_solution)
{
  Assert(
    gather_curvature,
    ExcMessage(
      "You are trying to reinit the curvature in a cell, but you did not enable curvature for the scratch data (gather_curvature=false)."));

  this->fe_values_curvature->reinit(curvature_cell);

  // Gather phase indicator gradient
  this->fe_values_curvature->get_function_values(current_curvature_solution,
                                                 this->curvature_values);
}



template <int dim>
void
NavierStokesScratchData<dim>::enable_void_fraction(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_void_fraction = true;
  // The hessian is needed to calculate the gradient of the velocity divergence
  // term in the stress tensor
  gather_hessian = true;

  // Contrary to the other physics, we enable the calculation of the JxW values
  // on the void fraction
  fe_values_void_fraction =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe,
                                    quadrature,
                                    update_values | update_gradients |
                                      update_JxW_values);
  point_evaluator_void_fraction =
    std::make_shared<FEPointEvaluation<1, dim>>(mapping, fe, update_values);

  // Void Fraction
  void_fraction_values = std::vector<double>(this->n_q_points);
  previous_void_fraction_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  void_fraction_gradient_values = std::vector<Tensor<1, dim>>(this->n_q_points);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_particle_field_projection(
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping,
  const FiniteElement<dim> &fe_particle_drag_proj,
  const FiniteElement<dim> &fe_particle_two_way_coupling_force_proj,
  const FiniteElement<dim> &fe_particle_velocity_proj,
  const FiniteElement<dim> &fe_particle_momentum_transfer_coefficient)
{
  gather_particle_field_project = true;

  fe_values_particle_drag = std::make_shared<FEValues<dim>>(
    mapping, fe_particle_drag_proj, quadrature, update_values);

  fe_values_particle_two_way_coupling_force =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_particle_two_way_coupling_force_proj,
                                    quadrature,
                                    update_values);
  fe_values_particle_velocity = std::make_shared<FEValues<dim>>(
    mapping, fe_particle_velocity_proj, quadrature, update_values);
  fe_values_particle_momentum_transfer_coefficient =
    std::make_shared<FEValues<dim>>(mapping,
                                    fe_particle_momentum_transfer_coefficient,
                                    quadrature,
                                    update_values);

  particle_drag_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_two_way_coupling_force_values =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  particle_momentum_transfer_coefficient_values =
    std::vector<double>(this->n_q_points);
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_void_fraction(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions)
{
  Assert(
    gather_void_fraction,
    ExcMessage(
      "You are trying to reinit the void fraction in a cell, but you did not enable void fraction for the scratch data (gather_void_fraction=false)."));

  this->fe_values_void_fraction->reinit(cell);

  // Gather void fraction (values, gradient)
  this->fe_values_void_fraction->get_function_values(
    current_solution, this->void_fraction_values);
  this->fe_values_void_fraction->get_function_gradients(
    current_solution, this->void_fraction_gradient_values);

  // Gather previous void fraction fraction values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values_void_fraction->get_function_values(
        previous_solutions[p], previous_void_fraction_values[p]);
    }
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_particle_fluid_interactions(
  const unsigned int n_global_max_particles_per_cell,
  const bool         enable_void_fraction_interpolation)
{
  gather_particles_information     = true;
  max_number_of_particles_per_cell = n_global_max_particles_per_cell;
  interpolated_void_fraction       = enable_void_fraction_interpolation;

  // Reinitialize vectors used to store flow information at the particle
  // location that do not rely on a quadrature
  particle_velocity =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  fluid_particle_relative_velocity_at_particle_location =
    std::vector<Tensor<1, dim>>(n_global_max_particles_per_cell);
  Re_particle = std::vector<double>(n_global_max_particles_per_cell);

  // This table is not used within an FeValues so we can avoid resizing it
  // dynamically
  density_at_particle_location.resize(n_global_max_particles_per_cell);
  kinematic_viscosity_at_particle_location.resize(
    n_global_max_particles_per_cell);
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_ale()
{
  gather_ale = true;
}

template <int dim>
void
NavierStokesScratchData<dim>::reinit_ale(const Parameters::ALE<dim> &ale)
{
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

      // Subtract u_ale from the advective velocity
      this->advective_velocity[q] -= velocity_ale;
    }
}

template <int dim>
void
NavierStokesScratchData<dim>::enable_mortar()
{
  gather_mortar = true;
}

template <int dim>
void
NavierStokesScratchData<dim>::reinit_mortar(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const Parameters::Mortar<dim>                        &mortar_parameters)
{
  // Get updated rotor angular velocity
  mortar_parameters.rotor_angular_velocity->set_time(
    this->simulation_control->get_current_time());
  const double rotor_angular_velocity =
    mortar_parameters.rotor_angular_velocity->value(Point<dim>());

  // Apply prescribed rotor angular velocity only at rotor cells
  // (material_id = 1)
  double omega = 0.0;
  if (cell->material_id() == 1)
    omega = rotor_angular_velocity;

  // The mortar implementation is not defined for more than two
  // materials. Throw when running debug mode as a security measure to prevent
  // this case
  Assert(
    cell->material_id() < 2,
    ExcMessage(
      "The material id in a cell was identified as equal or greater than 2, however the mortar implementation is not defined for more than two materials."));

  // Compute rotor linear velocity at quadrature points
  rotor_linear_velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      // Shift quadrature point by the center of rotation
      const auto p =
        fe_values.quadrature_point(q) - mortar_parameters.center_of_rotation;

      // Compute terms of u_ale
      rotor_linear_velocity_values[q] =
        LetheGridTools::angular_to_linear_velocity(
          omega, p, mortar_parameters.rotation_axis_direction);

      // Subtract u_ale from the advective velocity
      this->advective_velocity[q] -= rotor_linear_velocity_values[q];
    }
}


template <int dim>
void
NavierStokesScratchData<dim>::enable_heat_transfer(
  const FiniteElement<dim> &fe,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping)
{
  gather_temperature    = true;
  fe_values_temperature = std::make_shared<FEValues<dim>>(
    mapping, fe, quadrature, update_values | update_gradients);

  temperature_values = std::vector<double>(this->n_q_points);
  previous_temperature_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));
  temperature_gradients = std::vector<Tensor<1, dim>>(this->n_q_points);
  previous_temperature_gradients = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  fields.insert(
    std::pair<field, std::vector<double>>(field::temperature, n_q_points));
}


template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::reinit_heat_transfer(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions)
{
  this->fe_values_temperature->reinit(cell);

  // Gather temperature
  this->fe_values_temperature->get_function_values(current_solution,
                                                   this->temperature_values);

  // Gather temperature gradient
  this->fe_values_temperature->get_function_gradients(
    current_solution, this->temperature_gradients);

  // Extrapolate temperature and temperature gradient to t+dt using the BDF
  // if the simulation is transient
  const auto method = this->simulation_control->get_assembly_method();
  if (time_stepping_is_bdf(method))
    {
      // Gather previous temperature and temperature gradient values
      for (unsigned int p = 0; p < previous_solutions.size(); ++p)
        {
          fe_values_temperature->get_function_values(
            previous_solutions[p], this->previous_temperature_values[p]);
          fe_values_temperature->get_function_gradients(
            previous_solutions[p], this->previous_temperature_gradients[p]);
        }
      // Extrapolate temperature and temperature gradient
      std::vector<double> time_vector =
        this->simulation_control->get_simulation_times();
      bdf_extrapolate(time_vector,
                      previous_temperature_values,
                      number_of_previous_solutions(method),
                      temperature_values);
      bdf_extrapolate(time_vector,
                      previous_temperature_gradients,
                      number_of_previous_solutions(method),
                      temperature_gradients);
    }
}



template <int dim>
void
NavierStokesScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::temperature) &&
      gather_temperature)
    {
      set_field_vector(field::temperature,
                       this->temperature_values,
                       this->fields);
    }

  if (properties_manager.field_is_required(field::pressure))
    {
      set_field_vector(field::pressure, this->pressure_values, this->fields);
    }

  if (properties_manager.field_is_required(field::shear_rate))
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // Calculate shear rate (at each q)
          const Tensor<2, dim> shear_rate_tensor =
            velocity_gradients[q] + transpose(velocity_gradients[q]);

          // Calculate the shear rate magnitude
          shear_rate[q] = calculate_shear_rate_magnitude(shear_rate_tensor);
        }

      set_field_vector(field::shear_rate, shear_rate, this->fields);
    }

  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          // In the case of an incompressible flow, viscosity is the only
          // required property
          const auto rheology_model = properties_manager.get_rheology();
          rheology_model->vector_value(fields, kinematic_viscosity);

          rheology_model->get_kinematic_viscosity_for_stabilization_vector(
            fields, kinematic_viscosity_for_stabilization);

          kinematic_viscosity_scale =
            rheology_model->get_kinematic_viscosity_scale();

          if (!properties_manager.density_is_constant())
            {
              // For a weakly compressible flow, density variations will play a
              // role
              const auto density_model = properties_manager.get_density();
              density_model->vector_value(fields, density);
              density_psi = density_model->get_psi();
              density_ref = density_model->get_density_ref();
              // Dynamic viscosity is also necessary for compressible flows
              rheology_model->get_dynamic_viscosity_vector(density_ref,
                                                           fields,
                                                           dynamic_viscosity);
              rheology_model->get_dynamic_viscosity_for_stabilization_vector(
                density_ref, fields, dynamic_viscosity_for_stabilization);
            }
          else
            {
              density_scale = properties_manager.get_density_scale();
            }

          if (properties_manager.is_non_newtonian())
            {
              // Calculate derivative of viscosity with respect to shear rate
              rheology_model->vector_jacobian(
                fields, field::shear_rate, grad_kinematic_viscosity_shear_rate);
            }

          if (gather_temperature)
            {
              const auto thermal_expansion_model =
                properties_manager.get_thermal_expansion();
              thermal_expansion_model->vector_value(fields, thermal_expansion);
            }
          break;
        }
      case 2:
        {
          // We need both density and viscosity
          const auto density_model_0  = properties_manager.get_density(0);
          const auto rheology_model_0 = properties_manager.get_rheology(0);
          const auto density_model_1  = properties_manager.get_density(1);
          const auto rheology_model_1 = properties_manager.get_rheology(1);

          density_ref_0 = density_model_0->get_density_ref();
          density_ref_1 = density_model_1->get_density_ref();

          kinematic_viscosity_scale_0 =
            rheology_model_0->get_kinematic_viscosity_scale();
          kinematic_viscosity_scale_1 =
            rheology_model_1->get_kinematic_viscosity_scale();

          kinematic_viscosity_scale =
            std::max(kinematic_viscosity_scale_0, kinematic_viscosity_scale_1);

          // Gather properties from material interactions if necessary
          if (properties_manager.get_number_of_material_interactions() > 0)
            {
              const auto material_interaction_id =
                properties_manager.get_material_interaction_id(
                  material_interactions_type::fluid_fluid, 0, 1);
              // Gather surface tension
              const auto surface_tension_model =
                properties_manager.get_surface_tension(material_interaction_id);
              surface_tension_model->vector_value(fields, surface_tension);
              // Gather surface tension gradient only if necessary
              if (!properties_manager.surface_tension_is_constant())
                surface_tension_model->vector_jacobian(
                  fields, field::temperature, surface_tension_gradient);
            }

          density_model_0->vector_value(fields, density_0);
          rheology_model_0->get_dynamic_viscosity_vector(density_ref_0,
                                                         fields,
                                                         dynamic_viscosity_0);
          rheology_model_0->get_dynamic_viscosity_for_stabilization_vector(
            density_ref_0, fields, dynamic_viscosity_for_stabilization_0);

          density_model_1->vector_value(fields, density_1);
          rheology_model_1->get_dynamic_viscosity_vector(density_ref_1,
                                                         fields,
                                                         dynamic_viscosity_1);
          rheology_model_1->get_dynamic_viscosity_for_stabilization_vector(
            density_ref_1, fields, dynamic_viscosity_for_stabilization_1);
          if (gather_temperature)
            {
              const auto thermal_expansion_model_0 =
                properties_manager.get_thermal_expansion(0);
              const auto thermal_expansion_model_1 =
                properties_manager.get_thermal_expansion(1);
              thermal_expansion_model_0->vector_value(fields,
                                                      thermal_expansion_0);
              thermal_expansion_model_1->vector_value(fields,
                                                      thermal_expansion_1);
            }

          if (gather_cls && !gather_cahn_hilliard)
            {
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  double filtered_phase_value = this->filtered_phase_values[q];

                  density[q] = calculate_point_property(filtered_phase_value,
                                                        this->density_0[q],
                                                        this->density_1[q]);

                  dynamic_viscosity[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->dynamic_viscosity_0[q],
                                             this->dynamic_viscosity_1[q]);

                  kinematic_viscosity[q] = dynamic_viscosity[q] / density[q];

                  dynamic_viscosity_for_stabilization[q] =
                    calculate_point_property(
                      filtered_phase_value,
                      this->dynamic_viscosity_for_stabilization_0[q],
                      this->dynamic_viscosity_for_stabilization_1[q]);

                  thermal_expansion[q] =
                    calculate_point_property(filtered_phase_value,
                                             this->thermal_expansion_0[q],
                                             this->thermal_expansion_1[q]);

                  // Phase interpolate the product rho*beta per-phase to avoid
                  // cross term error from bilinear phase interpolation:
                  // interp_lin(rho)*interp_lin(beta) != interp_lin(rho*beta)
                  density_thermal_expansion[q] = calculate_point_property(
                    filtered_phase_value,
                    this->density_0[q] * this->thermal_expansion_0[q],
                    this->density_1[q] * this->thermal_expansion_1[q]);
                }

              // Gather density_psi for isothermal compressible NS equations
              if (!properties_manager.density_is_constant())
                {
                  density_psi_0 = density_model_0->get_psi();
                  density_psi_1 = density_model_1->get_psi();

                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      double filtered_phase_value =
                        this->filtered_phase_values[q];
                      // To avoid calculating twice in the assemblers
                      compressibility_multiplier[q] =
                        filtered_phase_value *
                          (density_0[q] / density_1[q] * density_psi_1 -
                           density_psi_0) +
                        density_psi_0;
                    }
                }

              for (unsigned p = 0; p < previous_phase_values.size(); ++p)
                {
                  // Blend the physical properties using the CLS field
                  for (unsigned int q = 0; q < this->n_q_points; ++q)
                    {
                      // Calculate previous density (right now assumes constant
                      // density model per phase)
                      previous_density[p][q] = calculate_point_property(
                        filter->filter_phase(this->previous_phase_values[p][q]),
                        this->density_0[q],
                        this->density_1[q]);
                    }
                }
              break;
            }
          else if (gather_cahn_hilliard && !gather_cls)
            {
              // Blend the physical properties using the CahnHilliard field
              for (unsigned int q = 0; q < this->n_q_points; ++q)
                {
                  double phase_order_cahn_hilliard_value =
                    this->filtered_phase_order_cahn_hilliard_values[q];

                  density[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->density_0[q],
                    this->density_1[q]);

                  dynamic_viscosity[q] = calculate_point_property_cahn_hilliard(
                    phase_order_cahn_hilliard_value,
                    this->dynamic_viscosity_0[q],
                    this->dynamic_viscosity_1[q]);

                  dynamic_viscosity_for_stabilization[q] =
                    calculate_point_property_cahn_hilliard(
                      phase_order_cahn_hilliard_value,
                      this->dynamic_viscosity_for_stabilization_0[q],
                      this->dynamic_viscosity_for_stabilization_1[q]);
                }
              break;
            }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}


template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::calculate_particle_fields_values(
  const typename DoFHandler<dim>::active_cell_iterator &particle_drag_cell,
  const typename DoFHandler<dim>::active_cell_iterator
    &particle_two_way_coupling_force_cell,
  const typename DoFHandler<dim>::active_cell_iterator &particle_velocity_cell,
  const typename DoFHandler<dim>::active_cell_iterator
                                 &particle_momentum_transfer_coefficient_cell,
  const VectorType               &particle_fluid_drag,
  const VectorType               &particle_fluid_force_two_way_coupling,
  const VectorType               &particle_velocity,
  const VectorType               &particle_momentum_transfer_coefficient,
  const Parameters::DragCoupling &drag_coupling)
{
  constexpr FEValuesExtractors::Vector vector_index(0);

  this->fe_values_particle_two_way_coupling_force->reinit(
    particle_two_way_coupling_force_cell);
  (*this->fe_values_particle_two_way_coupling_force)[vector_index]
    .get_function_values(particle_fluid_force_two_way_coupling,
                         this->particle_two_way_coupling_force_values);

  // particle_drag_values will remain zero in the implicit and semi-implicit
  // coupling, since the momentum transfer coefficient is used instead, while
  // particle_momentum_transfer_coefficient_values and
  // particle_velocity_values will remain zero in the fully explicit coupling
  if (drag_coupling == Parameters::DragCoupling::fully_explicit)
    {
      this->fe_values_particle_drag->reinit(particle_drag_cell);
      (*this->fe_values_particle_drag)[vector_index].get_function_values(
        particle_fluid_drag, this->particle_drag_values);
    }
  else
    {
      this->fe_values_particle_velocity->reinit(particle_velocity_cell);
      (*this->fe_values_particle_velocity)[vector_index].get_function_values(
        particle_velocity, this->particle_velocity_values);

      this->fe_values_particle_momentum_transfer_coefficient->reinit(
        particle_momentum_transfer_coefficient_cell);
      (*this->fe_values_particle_momentum_transfer_coefficient)
        .get_function_values(
          particle_momentum_transfer_coefficient,
          this->particle_momentum_transfer_coefficient_values);
    }
}


template <int dim>
template <typename PropertiesIndex>
void
NavierStokesScratchData<dim>::reinit_particle_fluid_forces()
{
  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();
      // Set the particle_fluid_interactions properties and vectors to 0
      for (int d = 0; d < dim; ++d)
        {
          particle_properties[PropertiesIndex::fem_force_two_way_coupling_x +
                              d]                                 = 0.;
          particle_properties[PropertiesIndex::fem_force_one_way_coupling_x +
                              d]                                 = 0.;
          particle_properties[PropertiesIndex::fem_drag_x + d]   = 0.;
          particle_properties[PropertiesIndex::fem_torque_x + d] = 0.;
          explicit_particle_volumetric_acceleration_on_fluid[d]  = 0.;
        }
    }
}

template <int dim>
template <typename PropertiesIndex>
void
NavierStokesScratchData<dim>::extract_particle_properties()
{
  this->average_particle_velocity = 0;
  // Loop over particles in cell
  this->total_particle_volume = 0;
  unsigned int i_particle     = 0;

  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();
      // Stores the values of particle velocity in a tensor
      particle_velocity[i_particle][0] =
        particle_properties[PropertiesIndex::v_x];
      particle_velocity[i_particle][1] =
        particle_properties[PropertiesIndex::v_y];
      if constexpr (dim == 3)
        particle_velocity[i_particle][2] =
          particle_properties[PropertiesIndex::v_z];

      if (!interpolated_void_fraction)
        total_particle_volume +=
          M_PI * pow(particle_properties[PropertiesIndex::dp], dim) / (2 * dim);

      average_particle_velocity += particle_velocity[i_particle];
      i_particle++;
    }
  number_of_particles = i_particle;
  if (number_of_particles != 0)
    { // Calculate the average particle velocity within the cell
      average_particle_velocity =
        average_particle_velocity / number_of_particles;
    }
}

template <int dim>
void
NavierStokesScratchData<dim>::calculate_cell_void_fraction(
  const double &total_particle_volume)
{
  cell_volume = compute_cell_measure_with_JxW(
    this->fe_values_void_fraction->get_JxW_values());

  // The cell_void_fraction vector needs to be resized to the correct size
  // since the quadrature used in the FEValues for the void fraction will
  // be tailored exactly to the number of particles.
  // We also initialize it to zero and it will be filled if a cell-averaged
  // void fraction is used.
  cell_void_fraction.resize(number_of_particles, 0);

  if (!this->interpolated_void_fraction)
    {
      double cell_void_fraction_bulk = 0;
      cell_void_fraction_bulk =
        (cell_volume - total_particle_volume) / cell_volume;

      for (unsigned int j = 0; j < number_of_particles; ++j)
        cell_void_fraction[j] = cell_void_fraction_bulk;
    }
}

template <int dim>
Quadrature<dim>
NavierStokesScratchData<dim>::gather_particles_reference_location()
{
  // Create local vector that will be used to spawn an in-situ quadrature to
  // interpolate at the locations of the particles
  std::vector<Point<dim>> particle_reference_location(number_of_particles);
  std::vector<double>     particle_weights(number_of_particles, 1);
  unsigned int            i_particle = 0;

  // Loop over particles in cell and cache their reference location
  for (auto &particle : pic)
    {
      // Store particle positions and weights
      // Reference location of the particle
      particle_reference_location[i_particle] =
        particle.get_reference_location();
      i_particle++;
    }

  // Return a quadrature for the Navier-Stokes equations that is based on the
  // particle reference location
  return Quadrature<dim>(particle_reference_location, particle_weights);
}


template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::calculate_fluid_fields_at_particle_location(
  const Quadrature<dim>                                &q_particles_location,
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  const VectorType               &present_velocity_pressure_solution,
  const VectorType               &previous_velocity_pressure_solution,
  const Parameters::DragCoupling &drag_coupling)
{
  FEValues<dim> fe_values_local_particles(this->fe_values.get_fe(),
                                          q_particles_location,
                                          update_gradients | update_values |
                                            update_hessians);

  // Reallocate memory for the fields to be interpolated at the particle
  // location This has to be done for every cell because deal.II expects the
  // size of the vector to be strictly equal to the number of points in the
  // quadrature
  fluid_velocity_at_particle_location.resize(number_of_particles);
  fluid_velocity_laplacian_at_particle_location.resize(number_of_particles);
  fluid_velocity_curls_at_particle_location_2d.resize(number_of_particles);
  fluid_velocity_curls_at_particle_location_3d.resize(number_of_particles);
  fluid_pressure_gradients_at_particle_location.resize(number_of_particles);

  // Take velocity_pressure_solution according to the type of coupling used.
  const auto &velocity_pressure_solution =
    drag_coupling == Parameters::DragCoupling::fully_implicit ?
      present_velocity_pressure_solution :
      previous_velocity_pressure_solution;

  fe_values_local_particles.reinit(velocity_cell);

  // Calculate all fluid properties at the particle location
  fe_values_local_particles[velocities].get_function_values(
    velocity_pressure_solution, fluid_velocity_at_particle_location);

  fe_values_local_particles[velocities].get_function_laplacians(
    velocity_pressure_solution, fluid_velocity_laplacian_at_particle_location);

  if constexpr (dim == 2)
    {
      fe_values_local_particles[velocities].get_function_curls(
        velocity_pressure_solution,
        fluid_velocity_curls_at_particle_location_2d);
    }
  else if constexpr (dim == 3)
    {
      fe_values_local_particles[velocities].get_function_curls(
        velocity_pressure_solution,
        fluid_velocity_curls_at_particle_location_3d);
    }

  fe_values_local_particles[pressure].get_function_gradients(
    velocity_pressure_solution, fluid_pressure_gradients_at_particle_location);
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::calculate_void_fraction_at_particle_location(
  const Quadrature<dim>                                &q_particles_location,
  const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
  const VectorType                                     &void_fraction_solution)
{
  Assert(
    gather_void_fraction,
    ExcMessage(
      "gather_void_fraction has been set to false, yet you are trying to gather the void fraction at the location of the particles. The scratch data is currently unaware of the finite element interpolation for the void fraction and the simulation will abort."));

  Assert(
    cell_void_fraction.size() == q_particles_location.size(),
    ExcMessage(
      "The vector for the void fraction at the particle location does not have the same size as the quadrature used to evaluate it."));

  Vector<double> local_void_fraction_dof_values(
    void_fraction_cell->get_fe().n_dofs_per_cell());
  void_fraction_cell->get_dof_values(void_fraction_solution,
                                     local_void_fraction_dof_values);

  point_evaluator_void_fraction->reinit(
    void_fraction_cell, make_array_view(q_particles_location.get_points()));
  point_evaluator_void_fraction->evaluate(
    make_array_view(local_void_fraction_dof_values), EvaluationFlags::values);

  for (unsigned int i_particle = 0; i_particle < q_particles_location.size();
       ++i_particle)
    cell_void_fraction[i_particle] =
      point_evaluator_void_fraction->get_value(i_particle);
}

template <int dim>
void
NavierStokesScratchData<dim>::calculate_fluid_properties_at_particle_location()
{
  AssertThrow(
    properties_manager.density_is_constant(),
    ExcMessage(
      "The lethe-fluid-particles solver only supports constant densities"));
  AssertThrow(
    !properties_manager.get_rheology()->is_non_newtonian_rheological_model() &&
      !properties_manager.get_rheology()->is_non_newtonian_rheological_model(),
    ExcMessage(
      "The lethe-fluid-particles solver only supports constant rheology model"));

  if (gather_cls)
    {
      for (unsigned int i_particle = 0; i_particle < number_of_particles;
           ++i_particle)
        {
          density_at_particle_location[i_particle] = calculate_point_property(
            filtered_phase_values_at_particle_location[i_particle],
            this->density_ref_0,
            this->density_ref_1);
          double dynamic_viscosity_at_particle_location =
            calculate_point_property(
              filtered_phase_values_at_particle_location[i_particle],
              this->kinematic_viscosity_scale_0 * this->density_ref_0,
              this->kinematic_viscosity_scale_1 * this->density_ref_1);

          kinematic_viscosity_at_particle_location[i_particle] =
            dynamic_viscosity_at_particle_location /
            density_at_particle_location[i_particle];
        }

      // Properties have been calculated, return
      return;
    }

  // Regular case without CLS
  for (unsigned int i_particle = 0; i_particle < number_of_particles;
       ++i_particle)
    {
      // Gather the kinematic viscosity and density at the particle
      // location assuming a constant kinematic viscosity and density in a
      // single fluid
      kinematic_viscosity_at_particle_location[i_particle] =
        kinematic_viscosity_scale;
      density_at_particle_location[i_particle] = density_scale;
    }
}

template <int dim>
template <typename PropertiesIndex>
void
NavierStokesScratchData<dim>::calculate_force_parameters_at_particle_location()
{
  // Tolerance for the calculation of the particle Reynolds number.
  // TODO -> Revisit if that is not a too high value.
  const double Re_particle_tolerance       = 1e-3;
  unsigned int i_particle                  = 0;
  average_fluid_particle_relative_velocity = 0;

  for (auto &particle : pic)
    {
      auto particle_properties = particle.get_properties();

      fluid_particle_relative_velocity_at_particle_location[i_particle] =
        fluid_velocity_at_particle_location[i_particle] -
        particle_velocity[i_particle];
      average_fluid_particle_relative_velocity +=
        fluid_particle_relative_velocity_at_particle_location[i_particle];

      Re_particle[i_particle] =
        Re_particle_tolerance +
        cell_void_fraction[i_particle] *
          fluid_particle_relative_velocity_at_particle_location[i_particle]
            .norm() *
          particle_properties[PropertiesIndex::dp] /
          (kinematic_viscosity_at_particle_location[i_particle] + DBL_MIN);
      i_particle++;
    }

  average_fluid_particle_relative_velocity =
    average_fluid_particle_relative_velocity / i_particle;
}

template <int dim>
template <typename VectorType>
void
NavierStokesScratchData<dim>::calculate_cls_at_particle_location(
  const Quadrature<dim>                                &q_particles_location,
  const typename DoFHandler<dim>::active_cell_iterator &phase_cell,
  const VectorType &current_filtered_solution)
{
  Assert(
    gather_cls,
    ExcMessage(
      "gather_cls has been set to false, yet you are trying to gather the CLS at the location of the particles. The scratch data is currently unaware of the finite element interpolation for CLS and the simulation will abort."));

  filtered_phase_values_at_particle_location.resize(number_of_particles);

  Vector<double> local_filtered_phase_dof_values(
    phase_cell->get_fe().n_dofs_per_cell());
  phase_cell->get_dof_values(current_filtered_solution,
                             local_filtered_phase_dof_values);

  point_evaluator_cls->reinit(
    phase_cell, make_array_view(q_particles_location.get_points()));
  point_evaluator_cls->evaluate(
    make_array_view(local_filtered_phase_dof_values), EvaluationFlags::values);

  for (unsigned int i_particle = 0; i_particle < number_of_particles;
       ++i_particle)
    filtered_phase_values_at_particle_location[i_particle] =
      point_evaluator_cls->get_value(i_particle);
}

template <int dim>
template <typename PropertiesIndex, typename VectorType>
void
NavierStokesScratchData<dim>::reinit_particle_fluid_interactions(
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
  const VectorType                      &present_velocity_pressure_solution,
  const VectorType                      &previous_velocity_pressure_solution,
  const VectorType                      &void_fraction_solution,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::DragCoupling        &drag_coupling)
{
  Assert(
    gather_particles_information,
    ExcMessage(
      "You are trying to reinit the fluid information at the particle location within a cell,"
      " but you did not enable the gathering of the information at the particle location for the scratch data (gather_particles_information=false)."));

  pic = particle_handler.particles_in_cell(velocity_cell);

  extract_particle_properties<PropertiesIndex>();
  reinit_particle_fluid_forces<PropertiesIndex>();

  calculate_cell_void_fraction(total_particle_volume);

  if (number_of_particles == 0)
    return;

  Quadrature<dim> q_particles_location = gather_particles_reference_location();

  calculate_fluid_fields_at_particle_location(
    q_particles_location,
    velocity_cell,
    present_velocity_pressure_solution,
    previous_velocity_pressure_solution,
    drag_coupling);

  if (this->interpolated_void_fraction)
    {
      calculate_void_fraction_at_particle_location(q_particles_location,
                                                   void_fraction_cell,
                                                   void_fraction_solution);
    }
  calculate_fluid_properties_at_particle_location();
  calculate_force_parameters_at_particle_location<PropertiesIndex>();
}

template <int dim>
template <typename PropertiesIndex, typename VectorType>
void
NavierStokesScratchData<dim>::reinit_particle_fluid_interactions(
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<dim>::active_cell_iterator &void_fraction_cell,
  const typename DoFHandler<dim>::active_cell_iterator &phase_cell,
  const VectorType                      &present_velocity_pressure_solution,
  const VectorType                      &previous_velocity_pressure_solution,
  const VectorType                      &void_fraction_solution,
  const Particles::ParticleHandler<dim> &particle_handler,
  const Parameters::DragCoupling        &drag_coupling,
  const VectorType                      &current_filtered_CLS_solution)
{
  pic = particle_handler.particles_in_cell(velocity_cell);

  extract_particle_properties<PropertiesIndex>();
  reinit_particle_fluid_forces<PropertiesIndex>();
  calculate_cell_void_fraction(total_particle_volume);

  if (number_of_particles == 0)
    return;

  Quadrature<dim> q_particles_location = gather_particles_reference_location();
  calculate_fluid_fields_at_particle_location(
    q_particles_location,
    velocity_cell,
    present_velocity_pressure_solution,
    previous_velocity_pressure_solution,
    drag_coupling);

  if (this->interpolated_void_fraction)
    {
      calculate_void_fraction_at_particle_location(q_particles_location,
                                                   void_fraction_cell,
                                                   void_fraction_solution);
    }
  calculate_cls_at_particle_location(q_particles_location,
                                     phase_cell,
                                     current_filtered_CLS_solution);
  calculate_fluid_properties_at_particle_location();
  calculate_force_parameters_at_particle_location<PropertiesIndex>();
}

template class NavierStokesScratchData<2>;
template class NavierStokesScratchData<3>;

template void
NavierStokesScratchData<2>::reinit_particle_fluid_forces<
  DEM::CFDDEMProperties::PropertiesIndex>();

template void
NavierStokesScratchData<2>::extract_particle_properties<
  DEM::CFDDEMProperties::PropertiesIndex>();

template void
NavierStokesScratchData<2>::reinit_particle_fluid_forces<
  DEM::CFDDEMMPProperties::PropertiesIndex>();

template void
NavierStokesScratchData<2>::extract_particle_properties<
  DEM::CFDDEMMPProperties::PropertiesIndex>();

template void
NavierStokesScratchData<3>::reinit_particle_fluid_forces<
  DEM::CFDDEMProperties::PropertiesIndex>();

template void
NavierStokesScratchData<3>::extract_particle_properties<
  DEM::CFDDEMProperties::PropertiesIndex>();

template void
NavierStokesScratchData<3>::reinit_particle_fluid_forces<
  DEM::CFDDEMMPProperties::PropertiesIndex>();

template void
NavierStokesScratchData<3>::extract_particle_properties<
  DEM::CFDDEMMPProperties::PropertiesIndex>();

template void
NavierStokesScratchData<2>::reinit<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions,
  const GlobalVectorType                             &sum_over_previous_stages,
  const std::shared_ptr<Function<2>>                 &forcing_function,
  const Tensor<1, 2>                                 &beta_force,
  const double                                        pressure_scaling_factor);

template void
NavierStokesScratchData<2>::reinit<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalBlockVectorType                        &current_solution,
  const std::vector<GlobalBlockVectorType>           &previous_solutions,
  const GlobalBlockVectorType                        &sum_over_previous_stages,
  const std::shared_ptr<Function<2>>                 &forcing_function,
  const Tensor<1, 2>                                 &beta_force,
  const double                                        pressure_scaling_factor);

template void
NavierStokesScratchData<2>::reinit_cls<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const GlobalVectorType                             &current_filtered_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<2>::reinit_projected_phase_indicator_gradient<
  GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator
                         &projected_phase_indicator_gradient_cell,
  const GlobalVectorType &current_projected_phase_indicator_gradient_solution);

template void
NavierStokesScratchData<2>::reinit_curvature<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &curvature_cell,
  const GlobalVectorType &current_curvature_solution);

template void
NavierStokesScratchData<2>::reinit_void_fraction<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<2>::reinit_heat_transfer<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<2>::reinit_cahn_hilliard<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const GlobalVectorType &current_filtered_solution);

template void
NavierStokesScratchData<2>::calculate_particle_fields_values<
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<2>::active_cell_iterator &particle_drag_cell,
  const typename DoFHandler<2>::active_cell_iterator
    &particle_two_way_coupling_force_cell,
  const typename DoFHandler<2>::active_cell_iterator &particle_velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator
    &particle_momentum_transfer_coefficient_cell,
  const LinearAlgebra::distributed::Vector<double> &particle_fluid_drag,
  const LinearAlgebra::distributed::Vector<double>
    &particle_fluid_force_two_way_coupling,
  const LinearAlgebra::distributed::Vector<double> &particle_velocity,
  const LinearAlgebra::distributed::Vector<double>
                                 &particle_momentum_transfer_coefficient,
  const Parameters::DragCoupling &drag_coupling);

template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling);

template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const typename DoFHandler<2>::active_cell_iterator &phase_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling,
  const GlobalVectorType              &current_filtered_CLS_solution);

// Called with deal.II vectors by the matrix-free CFD-DEM solvers through
// ParticleProjector::calculate_particle_fluid_forces_projection. With
// LETHE_USE_LDV, GlobalVectorType is already this type.
#ifndef LETHE_USE_LDV
template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const LinearAlgebra::distributed::Vector<double>
    &present_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double>
    &previous_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution,
  const Particles::ParticleHandler<2>              &particle_handler,
  const Parameters::DragCoupling                   &drag_coupling);
#endif

template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling);

template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const typename DoFHandler<2>::active_cell_iterator &phase_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<2> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling,
  const GlobalVectorType              &current_filtered_CLS_solution);

#ifndef LETHE_USE_LDV
template void
NavierStokesScratchData<2>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<2>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<2>::active_cell_iterator &void_fraction_cell,
  const LinearAlgebra::distributed::Vector<double>
    &present_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double>
    &previous_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution,
  const Particles::ParticleHandler<2>              &particle_handler,
  const Parameters::DragCoupling                   &drag_coupling);
#endif

template void
NavierStokesScratchData<3>::reinit<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions,
  const GlobalVectorType                             &sum_over_previous_stages,
  const std::shared_ptr<Function<3>>                 &forcing_function,
  const Tensor<1, 3>                                 &beta_force,
  const double                                        pressure_scaling_factor);

template void
NavierStokesScratchData<3>::reinit<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalBlockVectorType                        &current_solution,
  const std::vector<GlobalBlockVectorType>           &previous_solutions,
  const GlobalBlockVectorType                        &sum_over_previous_stages,
  const std::shared_ptr<Function<3>>                 &forcing_function,
  const Tensor<1, 3>                                 &beta_force,
  const double                                        pressure_scaling_factor);

template void
NavierStokesScratchData<3>::reinit_cls<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const GlobalVectorType                             &current_filtered_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<3>::reinit_projected_phase_indicator_gradient<
  GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator
                         &projected_phase_indicator_gradient_cell,
  const GlobalVectorType &current_projected_phase_indicator_gradient_solution);

template void
NavierStokesScratchData<3>::reinit_curvature<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &curvature_cell,
  const GlobalVectorType &current_curvature_solution);

template void
NavierStokesScratchData<3>::reinit_void_fraction<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<3>::reinit_heat_transfer<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions);

template void
NavierStokesScratchData<3>::reinit_cahn_hilliard<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const GlobalVectorType &current_filtered_solution);

template void
NavierStokesScratchData<3>::calculate_particle_fields_values<
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<3>::active_cell_iterator &particle_drag_cell,
  const typename DoFHandler<3>::active_cell_iterator
    &particle_two_way_coupling_force_cell,
  const typename DoFHandler<3>::active_cell_iterator &particle_velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator
    &particle_momentum_transfer_coefficient_cell,
  const LinearAlgebra::distributed::Vector<double> &particle_fluid_drag,
  const LinearAlgebra::distributed::Vector<double>
    &particle_fluid_force_two_way_coupling,
  const LinearAlgebra::distributed::Vector<double> &particle_velocity,
  const LinearAlgebra::distributed::Vector<double>
                                 &particle_momentum_transfer_coefficient,
  const Parameters::DragCoupling &drag_coupling);

template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling);

template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const typename DoFHandler<3>::active_cell_iterator &phase_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling,
  const GlobalVectorType              &current_filtered_CLS_solution);

#ifndef LETHE_USE_LDV
template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMProperties::PropertiesIndex,
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const LinearAlgebra::distributed::Vector<double>
    &present_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double>
    &previous_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution,
  const Particles::ParticleHandler<3>              &particle_handler,
  const Parameters::DragCoupling                   &drag_coupling);
#endif

template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling);

template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const typename DoFHandler<3>::active_cell_iterator &phase_cell,
  const GlobalVectorType              &present_velocity_pressure_solution,
  const GlobalVectorType              &previous_velocity_pressure_solution,
  const GlobalVectorType              &void_fraction_solution,
  const Particles::ParticleHandler<3> &particle_handler,
  const Parameters::DragCoupling      &drag_coupling,
  const GlobalVectorType              &current_filtered_CLS_solution);

#ifndef LETHE_USE_LDV
template void
NavierStokesScratchData<3>::reinit_particle_fluid_interactions<
  DEM::CFDDEMMPProperties::PropertiesIndex,
  LinearAlgebra::distributed::Vector<double>>(
  const typename DoFHandler<3>::active_cell_iterator &velocity_cell,
  const typename DoFHandler<3>::active_cell_iterator &void_fraction_cell,
  const LinearAlgebra::distributed::Vector<double>
    &present_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double>
    &previous_velocity_pressure_solution,
  const LinearAlgebra::distributed::Vector<double> &void_fraction_solution,
  const Particles::ParticleHandler<3>              &particle_handler,
  const Parameters::DragCoupling                   &drag_coupling);
#endif
