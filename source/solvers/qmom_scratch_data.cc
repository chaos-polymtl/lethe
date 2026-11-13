// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <core/bdf.h>
#include <core/vector.h>

#include <solvers/qmom_scratch_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <memory>
#include <string>
#include <vector>

template <int dim>
QMOMScratchData<dim>::QMOMScratchData(
  const FiniteElement<dim>  &fe_qmom,
  const Quadrature<dim>     &quadrature,
  const Quadrature<dim - 1> &face_quadrature,
  const Mapping<dim>        &mapping,
  const FiniteElement<dim>  &fe_fluid,
  const unsigned int         n_components)
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
  , velocities(0)
  , n_components(n_components)
{
  AssertThrow(fe_qmom.n_components() == 1,
              ExcMessage("The QMOM components share a scalar element, but the "
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
                         sd.fe_values_fluid.get_fe(),
                         sd.n_components)
{}

template <int dim>
void
QMOMScratchData<dim>::allocate()
{
  // Initialize size of arrays
  n_dofs           = fe_values.get_fe().n_dofs_per_cell();
  n_q_points       = fe_values.get_quadrature().size();
  n_face_q_points  = fe_interface_values.get_quadrature().size();
  n_interface_dofs = 0;

  // Initialize arrays related to quadrature
  JxW.resize(n_q_points);
  quadrature_points.resize(n_q_points);

  // Shape functions
  phi.reinit(n_q_points, n_dofs);
  grad_phi.reinit(n_q_points, n_dofs);

  // Previous solutions of the BDF schemes and source term of every component
  previous_values.assign(
    n_components,
    std::vector<std::vector<double>>(maximum_number_of_previous_solutions(),
                                     std::vector<double>(n_q_points)));
  source_values.assign(n_q_points, Vector<double>(n_components));

  // Velocity
  velocity_values.resize(n_q_points);

  // Arrays related to the faces. An interior face has at most the degrees of
  // freedom of its two cells.
  const unsigned int max_n_interface_dofs = 2 * n_dofs;
  face_phi_here.reinit(n_face_q_points, max_n_interface_dofs);
  face_phi_there.reinit(n_face_q_points, max_n_interface_dofs);
  face_jump_phi.reinit(n_face_q_points, max_n_interface_dofs);
  face_velocity_values.resize(n_face_q_points);
}

template <int dim>
void
QMOMScratchData<dim>::reinit(
  const typename DoFHandler<dim>::active_cell_iterator &cell)
{
  fe_values.reinit(cell);
  quadrature_points = fe_values.get_quadrature_points();

  for (unsigned int q = 0; q < n_q_points; ++q)
    {
      JxW[q] = fe_values.JxW(q);
      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          phi[q][i]      = fe_values.shape_value(i, q);
          grad_phi[q][i] = fe_values.shape_grad(i, q);
        }
    }
}

template <int dim>
void
QMOMScratchData<dim>::reinit_previous_solutions(
  const std::vector<std::vector<GlobalVectorType>> &previous_solutions)
{
  AssertDimension(previous_solutions.size(), n_components);

  for (unsigned int component = 0; component < n_components; ++component)
    for (unsigned int p = 0; p < previous_solutions[component].size(); ++p)
      fe_values.get_function_values(previous_solutions[component][p],
                                    previous_values[component][p]);
}

template <int dim>
void
QMOMScratchData<dim>::reinit_source(
  const QMOMFieldFunction<dim> &source_function)
{
  // All the components are evaluated at once at each quadrature point, since
  // the moments are computed together from the particle size distribution
  for (unsigned int q = 0; q < n_q_points; ++q)
    source_function.vector_value(quadrature_points[q], source_values[q]);
}

template <int dim>
template <typename VectorType>
void
QMOMScratchData<dim>::reinit_velocity(
  const typename DoFHandler<dim>::active_cell_iterator &fluid_cell,
  const VectorType                                     &fluid_solution)
{
  fe_values_fluid.reinit(fluid_cell);
  fe_values_fluid[velocities].get_function_values(fluid_solution,
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
  n_interface_dofs = fe_interface_values.n_current_interface_dofs();

  for (unsigned int q = 0; q < n_face_q_points; ++q)
    {
      // Gather the shape function quantities once per quadrature point, since
      // the FEInterfaceValues accessors map the interface dof to the cell dofs
      // at every call. The jump is the value on the "here" cell minus the
      // value on the "there" cell.
      for (unsigned int i = 0; i < n_interface_dofs; ++i)
        {
          face_phi_here[q][i]  = fe_interface_values.shape_value(true, i, q);
          face_phi_there[q][i] = fe_interface_values.shape_value(false, i, q);
          face_jump_phi[q][i]  = face_phi_here[q][i] - face_phi_there[q][i];
        }
    }
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
      fe_face_values_fluid[velocities].get_function_values(
        fluid_solution, face_velocity_values);
    }
  else
    {
      Assert(fe_subface_values_fluid != nullptr, ExcInternalError());
      fe_subface_values_fluid->reinit(fluid_cell, face_no, sub_face_no);
      (*fe_subface_values_fluid)[velocities].get_function_values(
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
