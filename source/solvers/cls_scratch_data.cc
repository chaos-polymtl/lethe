// SPDX-FileCopyrightText: Copyright (c) 2021-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/time_integration_utilities.h>
#include <core/utilities.h>
#include <core/vector.h>

#include <solvers/cls_scratch_data.h>

template <int dim>
CLSScratchData<dim>::CLSScratchData(
  const std::shared_ptr<SimulationControl> &simulation_control,
  const PhysicalPropertiesManager          &properties_manager,
  const FiniteElement<dim>                 &fe_cls,
  const Quadrature<dim>                    &quadrature,
  const Quadrature<dim - 1>                &face_quadrature,
  const Mapping<dim>                       &mapping,
  const FiniteElement<dim>                 &fe_fd)
  : simulation_control(simulation_control)
  , properties_manager(properties_manager)
  , fe_values_cls(mapping,
                  fe_cls,
                  quadrature,
                  update_values | update_gradients | update_quadrature_points |
                    update_hessians | update_JxW_values)
  , fe_interface_values_cls(mapping,
                            fe_cls,
                            face_quadrature,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_normal_vectors)
  , fe_values_fd(mapping, fe_fd, quadrature, update_values | update_gradients)
  , fe_face_values_fd(mapping, fe_fd, face_quadrature, update_values)
{
  allocate();
}

template <int dim>
CLSScratchData<dim>::CLSScratchData(const CLSScratchData<dim> &sd)
  : simulation_control(sd.simulation_control)
  , properties_manager(sd.properties_manager)
  , fe_values_cls(sd.fe_values_cls.get_mapping(),
                  sd.fe_values_cls.get_fe(),
                  sd.fe_values_cls.get_quadrature(),
                  update_values | update_gradients | update_quadrature_points |
                    update_hessians | update_JxW_values)
  , fe_interface_values_cls(sd.fe_interface_values_cls.get_mapping(),
                            sd.fe_interface_values_cls.get_fe(),
                            sd.fe_interface_values_cls.get_quadrature(),
                            update_values | update_quadrature_points |
                              update_JxW_values | update_normal_vectors)
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
CLSScratchData<dim>::allocate()
{
  // Initialize size of arrays
  this->n_q_points = fe_values_cls.get_quadrature().size();
  this->n_dofs     = fe_values_cls.get_fe().n_dofs_per_cell();

  // Initialize arrays related to quadrature
  this->JxW.reinit(n_q_points);

  // Initialize arrays related to velocity and pressure
  this->velocities_fd.first_vector_component = 0;
  // Velocity
  this->velocity_values = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_velocity_values = std::vector<std::vector<Tensor<1, dim>>>(
    maximum_number_of_previous_solutions(),
    std::vector<Tensor<1, dim>>(this->n_q_points));
  this->velocity_gradient_values =
    std::vector<Tensor<2, dim>>(this->n_q_points);
  this->velocity_divergences = std::vector<double>(n_q_points);

  this->present_phase_values = std::vector<double>(this->n_q_points);
  this->phase_gradients      = std::vector<Tensor<1, dim>>(this->n_q_points);
  this->previous_phase_gradients =
    std::vector<Tensor<1, dim>>(this->n_q_points);

  this->phase_laplacians = std::vector<double>(this->n_q_points);

  // Velocity for BDF schemes
  this->previous_phase_values =
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(this->n_q_points));

  // Initialize arrays related to shape functions
  this->phi.reinit(n_q_points, n_dofs);
  this->grad_phi.reinit(n_q_points, n_dofs);
  this->hess_phi.reinit(n_q_points, n_dofs);
  this->laplacian_phi.reinit(n_q_points, n_dofs);
}

template <int dim>
template <typename VectorType>
void
CLSScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions)
{
  fe_values_cls.reinit(cell);
  this->quadrature_points = fe_values_cls.get_quadrature_points();
  auto &fe_cls            = fe_values_cls.get_fe();

  // Compute cell diameter
  double cell_measure =
    compute_cell_measure_with_JxW(this->fe_values_cls.get_JxW_values());
  this->cell_size = compute_cell_diameter<dim>(cell_measure, fe_cls.degree);

  fe_values_cls.get_function_values(current_solution,
                                    this->present_phase_values);
  fe_values_cls.get_function_gradients(current_solution, this->phase_gradients);
  fe_values_cls.get_function_laplacians(current_solution,
                                        this->phase_laplacians);
  fe_values_cls.get_function_gradients(previous_solutions[0],
                                       this->previous_phase_gradients);


  // Gather previous cls values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      fe_values_cls.get_function_values(previous_solutions[p],
                                        this->previous_phase_values[p]);
    }

  for (unsigned int q = 0; q < this->n_q_points; ++q)
    {
      this->JxW[q] = fe_values_cls.JxW(q);

      for (unsigned int k = 0; k < this->n_dofs; ++k)
        {
          // Shape function
          this->phi[q][k]           = fe_values_cls.shape_value(k, q);
          this->grad_phi[q][k]      = fe_values_cls.shape_grad(k, q);
          this->hess_phi[q][k]      = fe_values_cls.shape_hessian(k, q);
          this->laplacian_phi[q][k] = trace(this->hess_phi[q][k]);
        }
    }
}

template <int dim>
template <typename VectorType>
void
CLSScratchData<dim>::reinit_internal_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                   &face_no,
  const unsigned int                                   &sub_face_no,
  const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
  const unsigned int                                   &neighbor_face_no,
  const unsigned int                                   &neighbor_sub_face_no,
  const VectorType                                     &current_solution)
{
  fe_interface_values_cls.reinit(cell,
                                 face_no,
                                 sub_face_no,
                                 neighbor_cell,
                                 neighbor_face_no,
                                 neighbor_sub_face_no);
  face_quadrature_points = fe_interface_values_cls.get_quadrature_points();

  n_interface_dofs = fe_interface_values_cls.n_current_interface_dofs();

  // BB TODO : Preallocate memory here
  values_here.resize(face_quadrature_points.size());
  values_there.resize(face_quadrature_points.size());
  phase_value_jump.resize(face_quadrature_points.size());

  fe_interface_values_cls.get_fe_face_values(0).get_function_values(
    current_solution, values_here);
  fe_interface_values_cls.get_fe_face_values(1).get_function_values(
    current_solution, values_there);

  fe_interface_values_cls.get_jump_in_function_values(current_solution,
                                                      phase_value_jump);
}

template <int dim>
template <typename VectorType>
void
CLSScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const VectorType                                     &current_solution,
  const std::vector<VectorType>                        &previous_solutions,
  const Parameters::ALE<dim>                           &ale)
{
  fe_values_fd.reinit(cell);

  fe_values_fd[velocities_fd].get_function_values(current_solution,
                                                  velocity_values);
  fe_values_fd[velocities_fd].get_function_gradients(current_solution,
                                                     velocity_gradient_values);

  for (unsigned int q = 0; q < this->n_q_points; ++q)
    {
      this->velocity_divergences[q] = trace(this->velocity_gradient_values[q]);
    }

  // Gather previous velocity values
  for (unsigned int p = 0; p < previous_solutions.size(); ++p)
    {
      fe_values_fd[velocities_fd].get_function_values(
        previous_solutions[p], this->previous_velocity_values[p]);
    }

  if (ale.enabled())
    {
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

          for (unsigned int p = 0; p < previous_solutions.size(); ++p)
            {
              this->previous_velocity_values[p][q] -= velocity_ale;
            }
        }
    }

  // Extrapolate velocity to t+dt using the BDF scheme if the simulation is
  // transient
  const auto method = this->simulation_control->get_assembly_method();
  if (time_stepping_is_bdf(method))
    {
      // Extrapolate velocity
      std::vector<double> time_vector =
        this->simulation_control->get_simulation_times();
      bdf_extrapolate(time_vector,
                      this->previous_velocity_values,
                      number_of_previous_solutions(method),
                      this->velocity_values);
    }
}

template <int dim>
template <typename VectorType>
void
CLSScratchData<dim>::reinit_face_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &velocity_cell,
  const unsigned int                                   &face_no,
  const VectorType                                     &velocity_solution,
  const Parameters::ALE<dim>                           &ale)
{
  fe_face_values_fd.reinit(velocity_cell, face_no);

  // BB note : Array could be pre-allocated
  face_velocity_values.resize(face_quadrature_points.size());

  fe_face_values_fd[velocities_fd].get_function_values(velocity_solution,
                                                       face_velocity_values);

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

template class CLSScratchData<2>;
template class CLSScratchData<3>;

template void
CLSScratchData<2>::reinit<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions);

template void
CLSScratchData<2>::reinit_internal_face<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const unsigned int                        &face_no,
  const unsigned int                        &sub_face_no,
  const DoFHandler<2>::active_cell_iterator &neighbor_cell,
  const unsigned int                        &neighbor_face_no,
  const unsigned int                        &neighbor_sub_face_no,
  const GlobalVectorType                    &current_solution);

template void
CLSScratchData<3>::reinit<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions);

template void
CLSScratchData<3>::reinit_internal_face<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const unsigned int                        &face_no,
  const unsigned int                        &sub_face_no,
  const DoFHandler<3>::active_cell_iterator &neighbor_cell,
  const unsigned int                        &neighbor_face_no,
  const unsigned int                        &neighbor_sub_face_no,
  const GlobalVectorType                    &current_solution);

template void
CLSScratchData<2>::reinit_velocity<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions,
  const Parameters::ALE<2>                  &ale);

template void
CLSScratchData<2>::reinit_face_velocity<GlobalVectorType>(
  const DoFHandler<2>::active_cell_iterator &velocity_cell,
  const unsigned int                        &face_no,
  const GlobalVectorType                    &velocity_solution,
  const Parameters::ALE<2>                  &ale);

template void
CLSScratchData<2>::reinit_velocity<GlobalBlockVectorType>(
  const DoFHandler<2>::active_cell_iterator &cell,
  const GlobalBlockVectorType               &current_solution,
  const std::vector<GlobalBlockVectorType>  &previous_solutions,
  const Parameters::ALE<2>                  &ale);

template void
CLSScratchData<2>::reinit_face_velocity<GlobalBlockVectorType>(
  const DoFHandler<2>::active_cell_iterator &velocity_cell,
  const unsigned int                        &face_no,
  const GlobalBlockVectorType               &velocity_solution,
  const Parameters::ALE<2>                  &ale);

template void
CLSScratchData<3>::reinit_velocity<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalVectorType                    &current_solution,
  const std::vector<GlobalVectorType>       &previous_solutions,
  const Parameters::ALE<3>                  &ale);

template void
CLSScratchData<3>::reinit_face_velocity<GlobalVectorType>(
  const DoFHandler<3>::active_cell_iterator &velocity_cell,
  const unsigned int                        &face_no,
  const GlobalVectorType                    &velocity_solution,
  const Parameters::ALE<3>                  &ale);

template void
CLSScratchData<3>::reinit_velocity<GlobalBlockVectorType>(
  const DoFHandler<3>::active_cell_iterator &cell,
  const GlobalBlockVectorType               &current_solution,
  const std::vector<GlobalBlockVectorType>  &previous_solutions,
  const Parameters::ALE<3>                  &ale);

template void
CLSScratchData<3>::reinit_face_velocity<GlobalBlockVectorType>(
  const DoFHandler<3>::active_cell_iterator &velocity_cell,
  const unsigned int                        &face_no,
  const GlobalBlockVectorType               &velocity_solution,
  const Parameters::ALE<3>                  &ale);
