// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/vector.h>

#include <solvers/tracer_scratch_data.h>

template <int dim>
TracerScratchData<dim>::TracerScratchData(
  const PhysicalPropertiesManager &properties_manager,
  const FiniteElement<dim>        &fe_tracer,
  const Quadrature<dim>           &quadrature,
  const Quadrature<dim - 1>       &face_quadrature,
  const Mapping<dim>              &mapping,
  const FiniteElement<dim>        &fe_fd,
  const bool                       uses_dg)
  : uses_dg(uses_dg)
  , properties_manager(properties_manager)
  , fe_values_tracer(mapping,
                     fe_tracer,
                     quadrature,
                     update_values | update_quadrature_points |
                       update_JxW_values | update_gradients |
                       (uses_dg ? update_default : update_hessians))
  , fe_interface_values_tracer(mapping,
                               fe_tracer,
                               face_quadrature,
                               update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values |
                                 update_normal_vectors)
  , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
  , fe_face_values_fd(mapping, fe_fd, face_quadrature, update_values)
{
  allocate();
}

template <int dim>
TracerScratchData<dim>::TracerScratchData(const TracerScratchData<dim> &sd)
  : uses_dg(sd.uses_dg)
  , properties_manager(sd.properties_manager)
  , fe_values_tracer(sd.fe_values_tracer.get_mapping(),
                     sd.fe_values_tracer.get_fe(),
                     sd.fe_values_tracer.get_quadrature(),
                     update_values | update_quadrature_points |
                       update_JxW_values | update_gradients |
                       (sd.uses_dg ? update_default : update_hessians))
  , fe_interface_values_tracer(sd.fe_interface_values_tracer.get_mapping(),
                               sd.fe_interface_values_tracer.get_fe(),
                               sd.fe_interface_values_tracer.get_quadrature(),
                               update_values | update_gradients |
                                 update_quadrature_points | update_JxW_values |
                                 update_normal_vectors)
  , fe_values_fd(sd.fe_values_fd.get_mapping(),
                 sd.fe_values_fd.get_fe(),
                 sd.fe_values_fd.get_quadrature(),
                 update_values | update_gradients)
  , fe_face_values_fd(sd.fe_face_values_fd.get_mapping(),
                      sd.fe_face_values_fd.get_fe(),
                      sd.fe_face_values_fd.get_quadrature(),
                      update_values)
{
  allocate();
}

template <int dim>
void
TracerScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_tracer.get_quadrature().size();
  this->n_dofs     = fe_values_tracer.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Forcing term array
  this->source = std::vector<double>(n_q_points);

  // Initialize arrays related to velocity
  this->velocities.first_vector_component = 0;
  // Velocity
  this->velocity_values      = std::vector<Tensor<1, dim>>(n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);

  // Tracer
  this->tracer_values    = std::vector<double>(n_q_points);
  this->tracer_gradients = std::vector<Tensor<1, dim>>(n_q_points);
  this->previous_tracer_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);
  this->tracer_laplacians              = std::vector<double>(n_q_points);
  this->tracer_diffusivity             = std::vector<double>(n_q_points);
  this->tracer_diffusivity_0           = std::vector<double>(n_q_points);
  this->tracer_diffusivity_1           = std::vector<double>(n_q_points);
  this->tracer_reaction_prefactor      = std::vector<double>(n_q_points);
  this->grad_tracer_reaction_prefactor = std::vector<double>(n_q_points);

  // Arrays related to the faces, which are only used by the DG assemblers
  const unsigned int n_face_q_points =
    fe_interface_values_tracer.get_quadrature().size();
  const unsigned int max_n_interface_dofs = 2 * n_dofs;
  this->tracer_diffusivity_face = std::vector<double>(n_face_q_points);
  this->face_velocity_values    = std::vector<Tensor<1, dim>>(n_face_q_points);
  this->values_here             = std::vector<double>(n_face_q_points);
  this->values_there            = std::vector<double>(n_face_q_points);
  this->tracer_value_jump       = std::vector<double>(n_face_q_points);
  this->gradients_here          = std::vector<Tensor<1, dim>>(n_face_q_points);
  this->tracer_average_gradient = std::vector<Tensor<1, dim>>(n_face_q_points);
  this->face_phi_here.reinit(n_face_q_points, max_n_interface_dofs);
  this->face_phi_there.reinit(n_face_q_points, max_n_interface_dofs);
  this->face_jump_phi.reinit(n_face_q_points, max_n_interface_dofs);
  this->face_average_grad_phi_dot_n.reinit(n_face_q_points,
                                           max_n_interface_dofs);

  // Solid signed distance function
  if (properties_manager.field_is_required(field::levelset))
    this->sdf_values = std::vector<double>(n_q_points);

  // Velocity for BDF schemes
  this->previous_tracer_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points));

  // Initialize arrays related to shape functions
  // Velocity shape functions
  this->phi.reinit(n_q_points, n_dofs);
  this->grad_phi.reinit(n_q_points, n_dofs);
  if (!uses_dg)
    {
      this->hess_phi.reinit(n_q_points, n_dofs);
      this->laplacian_phi.reinit(n_q_points, n_dofs);
    }

  if (properties_manager.field_is_required(field::levelset))
    fields.insert(
      std::pair<field, std::vector<double>>(field::levelset, n_q_points));
  if (properties_manager.field_is_required(field::tracer_concentration))
    fields.insert(
      std::pair<field, std::vector<double>>(field::tracer_concentration,
                                            n_q_points));
}


template <int dim>
template <typename VectorType>
void
TracerScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions,
  const Function<dim>                                  *source_function)
{
  this->fe_values_tracer.reinit(cell);

  quadrature_points = this->fe_values_tracer.get_quadrature_points();
  auto &fe_tracer   = this->fe_values_tracer.get_fe();

  source_function->value_list(quadrature_points, source);

  // Compute cell diameter
  double cell_measure =
    compute_cell_measure_with_JxW(this->fe_values_tracer.get_JxW_values());
  this->cell_size = compute_cell_diameter<dim>(cell_measure, fe_tracer.degree);

  // Gather tracer (values, gradient and laplacian)
  this->fe_values_tracer.get_function_values(current_solution,
                                             this->tracer_values);
  this->fe_values_tracer.get_function_gradients(current_solution,
                                                this->tracer_gradients);

  // The laplacian and the previous gradient are only used by the
  // stabilization of the continuous formulation
  if (!uses_dg)
    {
      this->fe_values_tracer.get_function_laplacians(current_solution,
                                                     this->tracer_laplacians);
      this->fe_values_tracer.get_function_gradients(
        previous_solutions[0], this->previous_tracer_gradients);
    }

  // Gather previous tracer values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      this->fe_values_tracer.get_function_values(previous_solutions[p],
                                                 previous_tracer_values[p]);
    }

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      this->JxW[q] = this->fe_values_tracer.JxW(q);

      for (unsigned int k = 0; k < n_dofs; ++k)
        {
          // Shape function
          this->phi[q][k]      = this->fe_values_tracer.shape_value(k, q);
          this->grad_phi[q][k] = this->fe_values_tracer.shape_grad(k, q);
        }

      if (!uses_dg)
        for (unsigned int k = 0; k < n_dofs; ++k)
          {
            this->hess_phi[q][k] = this->fe_values_tracer.shape_hessian(k, q);
            this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
          }
    }

  boundary_index = 0;
}


template <int dim>
void
TracerScratchData<dim>::reinit_signed_distance(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const std::shared_ptr<Shape<dim>>                    &immersed_solid_shape)
{
  if (properties_manager.field_is_required(field::levelset))
    {
      this->fe_values_tracer.reinit(cell);

      quadrature_points = this->fe_values_tracer.get_quadrature_points();

      Assert(
        immersed_solid_shape != nullptr,
        ExcMessage(
          "Shape is required for tracer assembly, but the shape is a nullptr"));

      for (unsigned int q = 0; q < n_q_points; q++)
        sdf_values[q] =
          immersed_solid_shape->value_with_cell_guess(quadrature_points[q],
                                                      cell);
    }
}


template <int dim>
void
TracerScratchData<dim>::reinit_signed_distance_at_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                   &face_no,
  const std::shared_ptr<Shape<dim>>                    &immersed_solid_shape)
{
  if (properties_manager.field_is_required(field::levelset))
    {
      fe_interface_values_tracer.reinit(cell, face_no);
      face_quadrature_points =
        fe_interface_values_tracer.get_quadrature_points();

      Assert(
        immersed_solid_shape != nullptr,
        ExcMessage(
          "Shape is required for tracer assembly, but the shape is a nullptr"));

      unsigned int n_face_q_points = face_quadrature_points.size();
      for (unsigned int q = 0; q < n_face_q_points; q++)
        sdf_values[q] =
          immersed_solid_shape->value_with_cell_guess(face_quadrature_points[q],
                                                      cell);
    }
}


template <int dim>
template <typename VectorType>
void
TracerScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator  &cell,
  const VectorType                                      &velocity_solution,
  const Parameters::ALE<dim>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<dim>> &drift_velocity)
{
  this->fe_values_fd.reinit(cell);

  this->fe_values_fd[velocities].get_function_values(velocity_solution,
                                                     velocity_values);

  this->fe_values_fd[velocities].get_function_divergences(velocity_solution,
                                                          velocity_divergences);

  // Add the drift velocity to the velocity to account for tracer drift flux
  // modeling
  Tensor<1, dim> drift_velocity_tensor;
  Vector<double> drift_velocity_vector(dim);

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      drift_velocity->vector_value(quadrature_points[q], drift_velocity_vector);
      for (int d = 0; d < dim; ++d)
        drift_velocity_tensor[d] = drift_velocity_vector[d];

      velocity_values[q] += drift_velocity_tensor;
    }

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
template <typename VectorType>
void
TracerScratchData<dim>::reinit_internal_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                   &face_no,
  const unsigned int                                   &sub_face_no,
  const typename DoFHandler<dim>::active_cell_iterator &neigh_cell,
  const unsigned int                                   &neigh_face_no,
  const unsigned int                                   &neigh_sub_face_no,
  const VectorType                                     &current_solution)
{
  fe_interface_values_tracer.reinit(
    cell, face_no, sub_face_no, neigh_cell, neigh_face_no, neigh_sub_face_no);
  face_quadrature_points = fe_interface_values_tracer.get_quadrature_points();

  n_interface_dofs = fe_interface_values_tracer.n_current_interface_dofs();

  const double extent_here = cell->measure() / cell->face(face_no)->measure();
  const double extent_there =
    neigh_cell->measure() / neigh_cell->face(neigh_face_no)->measure();

  penalty_factor = get_penalty_factor(fe_values_tracer.get_fe().degree,
                                      extent_here,
                                      extent_there);


  fe_interface_values_tracer.get_fe_face_values(0).get_function_values(
    current_solution, values_here);
  fe_interface_values_tracer.get_fe_face_values(1).get_function_values(
    current_solution, values_there);

  fe_interface_values_tracer.get_average_of_function_gradients(
    current_solution, tracer_average_gradient);

  const std::vector<Tensor<1, dim>> &normals =
    fe_interface_values_tracer.get_normal_vectors();

  for (unsigned int q = 0; q < face_quadrature_points.size(); ++q)
    {
      // The jump is defined as the value on cell 0 minus the value on cell 1
      tracer_value_jump[q] = values_here[q] - values_there[q];

      // Gather the shape function quantities once per quadrature point
      // since the FEInterfaceValues accessors map the interface dof to the
      // cell dofs at every call.
      for (unsigned int i = 0; i < n_interface_dofs; ++i)
        {
          face_phi_here[q][i] =
            fe_interface_values_tracer.shape_value(true, i, q);
          face_phi_there[q][i] =
            fe_interface_values_tracer.shape_value(false, i, q);
          face_jump_phi[q][i] = face_phi_here[q][i] - face_phi_there[q][i];
          face_average_grad_phi_dot_n[q][i] =
            fe_interface_values_tracer.average_of_shape_gradients(i, q) *
            normals[q];
        }
    }
}


template <int dim>
template <typename VectorType>
void
TracerScratchData<dim>::reinit_boundary_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                   &face_no,
  const unsigned int                                   &boundary_index,
  const VectorType                                     &current_solution)
{
  fe_interface_values_tracer.reinit(cell, face_no);
  const FEFaceValuesBase<dim> &fe_face =
    fe_interface_values_tracer.get_fe_face_values(0);
  face_quadrature_points = fe_interface_values_tracer.get_quadrature_points();

  fe_face.get_function_values(current_solution, values_here);
  fe_face.get_function_gradients(current_solution, gradients_here);

  this->boundary_index = boundary_index;

  const double extent_here = cell->measure() / cell->face(face_no)->measure();
  this->penalty_factor = get_penalty_factor(fe_values_tracer.get_fe().degree,
                                            extent_here,
                                            extent_here);
}


template <int dim>
template <typename VectorType>
void
TracerScratchData<dim>::reinit_face_velocity(
  const typename DoFHandler<dim>::active_cell_iterator  &velocity_cell,
  const unsigned int                                    &face_no,
  const VectorType                                      &velocity_solution,
  const Parameters::ALE<dim>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<dim>> &drift_velocity)
{
  fe_face_values_fd.reinit(velocity_cell, face_no);

  fe_face_values_fd[velocities].get_function_values(velocity_solution,
                                                    face_velocity_values);


  // Add the drift velocity to the velocity to account for tracer drift flux
  // modeling
  Tensor<1, dim> drift_velocity_tensor;
  Vector<double> drift_velocity_vector(dim);

  for (unsigned int q = 0; q < face_quadrature_points.size(); ++q)
    {
      drift_velocity->vector_value(face_quadrature_points[q],
                                   drift_velocity_vector);
      for (int d = 0; d < dim; ++d)
        drift_velocity_tensor[d] = drift_velocity_vector[d];

      face_velocity_values[q] += drift_velocity_tensor;
    }


  if (!ale.enabled())
    return;

  // ALE enabled, so extract the ALE velocity and subtract it from the
  // velocity obtained from the fluid dynamics
  Tensor<1, dim>                                  velocity_ale;
  std::shared_ptr<Functions::ParsedFunction<dim>> velocity_ale_function =
    ale.velocity;
  Vector<double> velocity_ale_vector(dim);

  for (unsigned int q = 0; q < face_quadrature_points.size(); ++q)
    {
      velocity_ale_function->vector_value(face_quadrature_points[q],
                                          velocity_ale_vector);
      for (int d = 0; d < dim; ++d)
        velocity_ale[d] = velocity_ale_vector[d];

      face_velocity_values[q] -= velocity_ale;
    }
}


template <int dim>
void
TracerScratchData<dim>::reinit_mortar(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const Parameters::Mortar<dim>                        &mortar_parameters,
  const double                                         &time)
{
  // Get updated rotor angular velocity
  mortar_parameters.rotor_angular_velocity->set_time(time);
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
      const auto p = fe_values_tracer.quadrature_point(q) -
                     mortar_parameters.center_of_rotation;

      // Compute terms of u_ale
      rotor_linear_velocity_values[q] =
        LetheGridTools::angular_to_linear_velocity(
          omega, p, mortar_parameters.rotation_axis_direction);

      // Subtract u_ale from the advective velocity
      velocity_values[q] -= rotor_linear_velocity_values[q];
    }
}


template <int dim>
void
TracerScratchData<dim>::calculate_face_physical_properties()
{
  const auto diffusivity_model = properties_manager.get_tracer_diffusivity();

  // If the tracer diffusivity depends on the level set, take this into
  // account
  if (properties_manager.field_is_required(field::levelset))
    set_field_vector(field::levelset, sdf_values, fields);


  diffusivity_model->vector_value(fields, tracer_diffusivity_face);
}


template <int dim>
void
TracerScratchData<dim>::calculate_physical_properties()
{
  if (properties_manager.field_is_required(field::levelset))
    set_field_vector(field::levelset, this->sdf_values, this->fields);

  if (properties_manager.field_is_required(field::tracer_concentration))
    set_field_vector(field::tracer_concentration,
                     this->tracer_values,
                     this->fields);

  switch (properties_manager.get_number_of_fluids())
    {
      case 1:
        {
          const auto diffusivity_model =
            properties_manager.get_tracer_diffusivity();
          diffusivity_model->vector_value(fields, tracer_diffusivity);
          const auto reaction_prefactor_model =
            properties_manager.get_tracer_reaction_prefactor();
          reaction_prefactor_model->vector_value(fields,
                                                 tracer_reaction_prefactor);
          reaction_prefactor_model->vector_jacobian(
            fields,
            field::tracer_concentration,
            grad_tracer_reaction_prefactor);

          break;
        }
      case 2:
        {
          throw std::runtime_error(
            "Unsupported number of fluids (2) with the tracer physics.");
          // In this case, we need both density and viscosity
          const auto diffusivity_models =
            properties_manager.get_tracer_diffusivity_vector();

          diffusivity_models[0]->vector_value(fields, tracer_diffusivity_0);
          diffusivity_models[1]->vector_value(fields, tracer_diffusivity_1);

          // TODO Incomplete at the present time because the tracer CLS
          // complete is not finished Blend the physical properties
          // using the CLS field
          // for (unsigned int q = 0; q < this->n_q_points; ++q)
          // {
          //          tracer_diffusivity[q] =
          //            calculate_point_property(this->phase_values[q],
          //                                     this->density_0[q],
          //                                     this->density_1[q]);
          // }
          break;
        }
      default:
        throw std::runtime_error("Unsupported number of fluids (>2)");
    }
}


template class TracerScratchData<2>;
template class TracerScratchData<3>;

template void
TracerScratchData<2>::reinit<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions,
  const Function<2>                                  *source_function);

template void
TracerScratchData<2>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator  &cell,
  const GlobalVectorType                              &velocity_solution,
  const Parameters::ALE<2>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<2>> &drift_velocity);

template void
TracerScratchData<2>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator  &cell,
  const GlobalBlockVectorType                         &velocity_solution,
  const Parameters::ALE<2>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<2>> &drift_velocity);

template void
TracerScratchData<2>::reinit_internal_face<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const unsigned int                                 &face_no,
  const unsigned int                                 &sub_face_no,
  const typename DoFHandler<2>::active_cell_iterator &neigh_cell,
  const unsigned int                                 &neigh_face_no,
  const unsigned int                                 &neigh_sub_face_no,
  const GlobalVectorType                             &current_solution);

template void
TracerScratchData<2>::reinit_boundary_face<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &cell,
  const unsigned int                                 &face_no,
  const unsigned int                                 &boundary_index,
  const GlobalVectorType                             &current_solution);

template void
TracerScratchData<2>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator  &velocity_cell,
  const unsigned int                                  &face_no,
  const GlobalVectorType                              &velocity_solution,
  const Parameters::ALE<2>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<2>> &drift_velocity);

template void
TracerScratchData<2>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator  &velocity_cell,
  const unsigned int                                  &face_no,
  const GlobalBlockVectorType                         &velocity_solution,
  const Parameters::ALE<2>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<2>> &drift_velocity);

template void
TracerScratchData<3>::reinit<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                             &current_solution,
  const std::vector<GlobalVectorType>                &previous_solutions,
  const Function<3>                                  *source_function);

template void
TracerScratchData<3>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator  &cell,
  const GlobalVectorType                              &velocity_solution,
  const Parameters::ALE<3>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<3>> &drift_velocity);

template void
TracerScratchData<3>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator  &cell,
  const GlobalBlockVectorType                         &velocity_solution,
  const Parameters::ALE<3>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<3>> &drift_velocity);

template void
TracerScratchData<3>::reinit_internal_face<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const unsigned int                                 &face_no,
  const unsigned int                                 &sub_face_no,
  const typename DoFHandler<3>::active_cell_iterator &neigh_cell,
  const unsigned int                                 &neigh_face_no,
  const unsigned int                                 &neigh_sub_face_no,
  const GlobalVectorType                             &current_solution);

template void
TracerScratchData<3>::reinit_boundary_face<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &cell,
  const unsigned int                                 &face_no,
  const unsigned int                                 &boundary_index,
  const GlobalVectorType                             &current_solution);

template void
TracerScratchData<3>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator  &velocity_cell,
  const unsigned int                                  &face_no,
  const GlobalVectorType                              &velocity_solution,
  const Parameters::ALE<3>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<3>> &drift_velocity);

template void
TracerScratchData<3>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator  &velocity_cell,
  const unsigned int                                  &face_no,
  const GlobalBlockVectorType                         &velocity_solution,
  const Parameters::ALE<3>                            &ale,
  const std::shared_ptr<Functions::ParsedFunction<3>> &drift_velocity);
