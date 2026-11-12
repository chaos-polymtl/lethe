// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/vector.h>

#include <solvers/qmom_scratch_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <array>
#include <memory>

template <int dim>
QMOMScratchData<dim>::QMOMScratchData(
  const FiniteElement<dim>  &fe_qmom,
  const Quadrature<dim>     &quadrature,
  const Quadrature<dim - 1> &face_quadrature,
  const Mapping<dim>        &mapping,
  const FiniteElement<dim>  &fe_fluid)
  : fe_values(mapping,
              fe_qmom,
              quadrature,
              update_values | update_gradients | update_quadrature_points |
                update_JxW_values)
  , fe_interface_values(mapping,
                        fe_qmom,
                        face_quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values | update_normal_vectors)
  , fe_values_fluid(mapping, fe_fluid, quadrature, update_values)
  , fe_face_values_fluid(mapping, fe_fluid, face_quadrature, update_values)
  , velocity_extractor(0)
{
  // Hanging nodes, and therefore subfaces, only exist on hypercube meshes
  if (fe_qmom.reference_cell().is_hyper_cube())
    fe_subface_values_fluid = std::make_unique<FESubfaceValues<dim>>(
      mapping, fe_fluid, face_quadrature, update_values);

  allocate();
}

template <int dim>
QMOMScratchData<dim>::QMOMScratchData(const QMOMScratchData<dim> &sd)
  : QMOMScratchData<dim>(sd.fe_values.get_fe(),
                         sd.fe_values.get_quadrature(),
                         sd.fe_interface_values.get_quadrature(),
                         sd.fe_values.get_mapping(),
                         sd.fe_values_fluid.get_fe())
{}

template <int dim>
void
QMOMScratchData<dim>::allocate()
{
  const FiniteElement<dim> &fe = fe_values.get_fe();

  n_components    = fe.n_components();
  n_dofs          = fe.n_dofs_per_cell();
  n_q_points      = fe_values.get_quadrature().size();
  n_face_q_points = fe_interface_values.get_quadrature().size();

  // All the components of the QMOM FESystem are primitive, so that each
  // degree of freedom belongs to a single component
  component_of_dof.resize(n_dofs);
  for (unsigned int i = 0; i < n_dofs; ++i)
    component_of_dof[i] = fe.system_to_component_index(i).first;

  values.assign(n_q_points, Vector<double>(n_components));
  velocity_values.resize(n_q_points);
  previous_velocity_values.resize(n_q_points);

  n_interface_dofs = 0;
  component_of_interface_dof.reserve(2 * n_dofs);
  values_here.assign(n_face_q_points, Vector<double>(n_components));
  values_there.assign(n_face_q_points, Vector<double>(n_components));
  face_velocity_values.resize(n_face_q_points);
  previous_face_velocity_values.resize(n_face_q_points);
  upwind_values.reinit(n_components);
}

template <int dim>
void
QMOMScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const GlobalVectorType                               &solution)
{
  fe_values.reinit(cell);
  fe_values.get_function_values(solution, values);
}

template <int dim>
template <typename VectorType>
void
QMOMScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
  const VectorType                                     &previous_fluid_solution,
  const VectorType                                     &present_fluid_solution,
  const double                                          interpolation_factor)
{
  fe_values_fluid.reinit(fluid_cell);
  fe_values_fluid[velocity_extractor].get_function_values(
    previous_fluid_solution, previous_velocity_values);
  fe_values_fluid[velocity_extractor].get_function_values(
    present_fluid_solution, velocity_values);

  for (unsigned int q = 0; q < n_q_points; ++q)
    velocity_values[q] =
      (1. - interpolation_factor) * previous_velocity_values[q] +
      interpolation_factor * velocity_values[q];
}

template <int dim>
void
QMOMScratchData<dim>::reinit_internal_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                    face_no,
  const unsigned int                                    sub_face_no,
  const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
  const unsigned int                                    neighbor_face_no,
  const unsigned int                                    neighbor_sub_face_no,
  const GlobalVectorType                               &solution)
{
  fe_interface_values.reinit(cell,
                             face_no,
                             sub_face_no,
                             neighbor_cell,
                             neighbor_face_no,
                             neighbor_sub_face_no);

  fe_interface_values.get_fe_face_values(0).get_function_values(solution,
                                                                values_here);
  fe_interface_values.get_fe_face_values(1).get_function_values(solution,
                                                                values_there);

  // Both cells share the same FESystem, so the component of an interface
  // degree of freedom can be read from either side
  const FiniteElement<dim> &fe = fe_values.get_fe();
  n_interface_dofs             = fe_interface_values.n_current_interface_dofs();
  component_of_interface_dof.resize(n_interface_dofs);
  for (unsigned int i = 0; i < n_interface_dofs; ++i)
    {
      const std::array<unsigned int, 2> local_dofs =
        fe_interface_values.interface_dof_to_dof_indices(i);
      const unsigned int local_dof =
        local_dofs[0] != numbers::invalid_unsigned_int ? local_dofs[0] :
                                                         local_dofs[1];
      component_of_interface_dof[i] =
        fe.system_to_component_index(local_dof).first;
    }
}

template <int dim>
void
QMOMScratchData<dim>::reinit_boundary_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                    face_no,
  const GlobalVectorType                               &solution)
{
  fe_interface_values.reinit(cell, face_no);
  fe_interface_values.get_fe_face_values(0).get_function_values(solution,
                                                                values_here);
}

template <int dim>
template <typename VectorType>
void
QMOMScratchData<dim>::reinit_face_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
  const unsigned int                                    face_no,
  const unsigned int                                    sub_face_no,
  const VectorType                                     &previous_fluid_solution,
  const VectorType                                     &present_fluid_solution,
  const double                                          interpolation_factor)
{
  // The fluid velocity is continuous: evaluating it from the first side of
  // the face, at the quadrature points of fe_interface_values, is sufficient.
  // On a refined face, these points lie on the subface.
  const FEFaceValuesBase<dim> *fluid_face_values = nullptr;
  if (sub_face_no == numbers::invalid_unsigned_int)
    {
      fe_face_values_fluid.reinit(fluid_cell, face_no);
      fluid_face_values = &fe_face_values_fluid;
    }
  else
    {
      Assert(fe_subface_values_fluid != nullptr, ExcInternalError());
      fe_subface_values_fluid->reinit(fluid_cell, face_no, sub_face_no);
      fluid_face_values = fe_subface_values_fluid.get();
    }

  (*fluid_face_values)[velocity_extractor].get_function_values(
    previous_fluid_solution, previous_face_velocity_values);
  (*fluid_face_values)[velocity_extractor].get_function_values(
    present_fluid_solution, face_velocity_values);

  for (unsigned int q = 0; q < n_face_q_points; ++q)
    face_velocity_values[q] =
      (1. - interpolation_factor) * previous_face_velocity_values[q] +
      interpolation_factor * face_velocity_values[q];
}

template class QMOMScratchData<2>;
template class QMOMScratchData<3>;

template void
QMOMScratchData<2>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const GlobalVectorType &,
  const GlobalVectorType &,
  const double);
template void
QMOMScratchData<3>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const GlobalVectorType &,
  const GlobalVectorType &,
  const double);
template void
QMOMScratchData<2>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const GlobalBlockVectorType &,
  const GlobalBlockVectorType &,
  const double);
template void
QMOMScratchData<3>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const GlobalBlockVectorType &,
  const GlobalBlockVectorType &,
  const double);

template void
QMOMScratchData<2>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalVectorType &,
  const GlobalVectorType &,
  const double);
template void
QMOMScratchData<3>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalVectorType &,
  const GlobalVectorType &,
  const double);
template void
QMOMScratchData<2>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalBlockVectorType &,
  const GlobalBlockVectorType &,
  const double);
template void
QMOMScratchData<3>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalBlockVectorType &,
  const GlobalBlockVectorType &,
  const double);
