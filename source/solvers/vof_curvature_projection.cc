// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_curvature_projection.h>

template <int dim>
void
VOFCurvatureProjection<dim>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFCurvatureProjectionScratchData<dim>             &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point);

  const DoFHandler<dim> *dof_handler_phase_gradient_projection =
    this->subequations->get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);

  typename DoFHandler<dim>::active_cell_iterator phase_gradient_projection_cell(
    &(*this->triangulation),
    cell->level(),
    cell->index(),
    dof_handler_phase_gradient_projection);

  scratch_data.reinit_projected_phase_gradient(
    phase_gradient_projection_cell,
    *this->subequations->get_solution(
      VOFSubequationsID::phase_gradient_projection));

  copy_data.reset();

  this->assembler->assemble_matrix(scratch_data, copy_data);

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim>
void
VOFCurvatureProjection<dim>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  VOFCurvatureProjectionScratchData<dim>                &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point);

  const DoFHandler<dim> *dof_handler_phase_gradient_projection =
    this->subequations->get_dof_handler(
      VOFSubequationsID::phase_gradient_projection);

  typename DoFHandler<dim>::active_cell_iterator phase_gradient_projection_cell(
    &(*this->triangulation),
    cell->level(),
    cell->index(),
    dof_handler_phase_gradient_projection);

  scratch_data.reinit_projected_phase_gradient(
    phase_gradient_projection_cell,
    *this->subequations->get_solution(
      VOFSubequationsID::phase_gradient_projection));

  copy_data.reset();

  this->assembler->assemble_rhs(scratch_data, copy_data);

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template class VOFCurvatureProjection<2>;
template class VOFCurvatureProjection<3>;
