// SPDX-FileCopyrightText: Copyright (c) 2020-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_phase_gradient_projection.h>

template <int dim, typename ScratchDataType>
void
VOFPhaseGradientProjection<dim, ScratchDataType>::assemble_local_system_matrix(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  ScratchDataType                                      &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point);

  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics->get_dof_handler(PhysicsID::VOF);

  typename DoFHandler<dim>::active_cell_iterator vof_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_vof);

  scratch_data.reinit_vof(
    vof_cell, *this->multiphysics->get_filtered_solution(PhysicsID::VOF));

  copy_data.reset();

  this->assembler->assemble_matrix(scratch_data, copy_data);

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template <int dim, typename ScratchDataType>
void
VOFPhaseGradientProjection<dim, ScratchDataType>::assemble_local_system_rhs(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  ScratchDataType                                      &scratch_data,
  StabilizedMethodsCopyData                            &copy_data)
{
  copy_data.cell_is_local = cell->is_locally_owned();
  if (!cell->is_locally_owned())
    return;

  scratch_data.reinit(cell, this->evaluation_point);

  const DoFHandler<dim> *dof_handler_vof =
    this->multiphysics->get_dof_handler(PhysicsID::VOF);

  typename DoFHandler<dim>::active_cell_iterator vof_cell(
    &(*this->triangulation), cell->level(), cell->index(), dof_handler_vof);

  scratch_data.reinit_vof(
    vof_cell, *this->multiphysics->get_filtered_solution(PhysicsID::VOF));

  copy_data.reset();

  this->assembler->assemble_rhs(scratch_data, copy_data);

  cell->get_dof_indices(copy_data.local_dof_indices);
}

template class VOFPhaseGradientProjection<
  2,
  VOFPhaseGradientProjectionScratchData<2>>;
template class VOFPhaseGradientProjection<
  3,
  VOFPhaseGradientProjectionScratchData<3>>;
