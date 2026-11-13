// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/vector.h>

#include <solvers/qmom_scratch_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <memory>
#include <string>

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
              update_values | update_gradients | update_JxW_values)
  , fe_interface_values(mapping,
                        fe_qmom,
                        face_quadrature,
                        update_values | update_quadrature_points |
                          update_JxW_values | update_normal_vectors)
  , fe_values_fluid(mapping, fe_fluid, quadrature, update_values)
  , fe_face_values_fluid(mapping, fe_fluid, face_quadrature, update_values)
  , velocity_extractor(0)
{
  AssertThrow(fe_qmom.n_components() == 1,
              ExcMessage("The QMOM fields share a scalar element, but the "
                         "element given to the scratch data has " +
                         std::to_string(fe_qmom.n_components()) +
                         " components."));

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
  n_dofs          = fe_values.get_fe().n_dofs_per_cell();
  n_q_points      = fe_values.get_quadrature().size();
  n_face_q_points = fe_interface_values.get_quadrature().size();

  JxW.resize(n_q_points);
  phi.reinit(n_dofs, n_q_points);
  grad_phi.reinit(n_dofs, n_q_points);
  velocity_values.resize(n_q_points);
  face_velocity_values.resize(n_face_q_points);
}

template <int dim>
void
QMOMScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  fe_values.reinit(cell);

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      JxW[q] = fe_values.JxW(q);
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          phi[i][q]      = fe_values.shape_value(i, q);
          grad_phi[i][q] = fe_values.shape_grad(i, q);
        }
    }
}

template <int dim>
template <typename VectorType>
void
QMOMScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
  const VectorType                                     &fluid_solution)
{
  fe_values_fluid.reinit(fluid_cell);
  fe_values_fluid[velocity_extractor].get_function_values(fluid_solution,
                                                          velocity_values);
}

template <int dim>
void
QMOMScratchData<dim>::reinit_internal_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                    face_no,
  const unsigned int                                    sub_face_no,
  const typename DoFHandler<dim>::active_cell_iterator &neighbor_cell,
  const unsigned int                                    neighbor_face_no,
  const unsigned int                                    neighbor_sub_face_no)
{
  fe_interface_values.reinit(cell,
                             face_no,
                             sub_face_no,
                             neighbor_cell,
                             neighbor_face_no,
                             neighbor_sub_face_no);
}

template <int dim>
void
QMOMScratchData<dim>::reinit_boundary_face(
  const typename DoFHandler<dim>::active_cell_iterator &cell,
  const unsigned int                                    face_no)
{
  fe_interface_values.reinit(cell, face_no);
}

template <int dim>
template <typename VectorType>
void
QMOMScratchData<dim>::reinit_face_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
  const unsigned int                                    face_no,
  const unsigned int                                    sub_face_no,
  const VectorType                                     &fluid_solution)
{
  // The fluid velocity is continuous: evaluating it from the first side of
  // the face, at the quadrature points of fe_interface_values, is sufficient.
  // On a refined face, these points lie on the subface.
  if (sub_face_no == numbers::invalid_unsigned_int)
    {
      fe_face_values_fluid.reinit(fluid_cell, face_no);
      fe_face_values_fluid[velocity_extractor].get_function_values(
        fluid_solution, face_velocity_values);
    }
  else
    {
      Assert(fe_subface_values_fluid != nullptr, ExcInternalError());
      fe_subface_values_fluid->reinit(fluid_cell, face_no, sub_face_no);
      (*fe_subface_values_fluid)[velocity_extractor].get_function_values(
        fluid_solution, face_velocity_values);
    }
}

template class QMOMScratchData<2>;
template class QMOMScratchData<3>;

template void
QMOMScratchData<2>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const GlobalVectorType &);
template void
QMOMScratchData<3>::reinit_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const GlobalVectorType &);
template void
QMOMScratchData<2>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const GlobalBlockVectorType &);
template void
QMOMScratchData<3>::reinit_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const GlobalBlockVectorType &);

template void
QMOMScratchData<2>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalVectorType &);
template void
QMOMScratchData<3>::reinit_face_velocity<GlobalVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalVectorType &);
template void
QMOMScratchData<2>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<2>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalBlockVectorType &);
template void
QMOMScratchData<3>::reinit_face_velocity<GlobalBlockVectorType>(
  const typename DoFHandler<3>::active_cell_iterator &,
  const unsigned int,
  const unsigned int,
  const GlobalBlockVectorType &);
