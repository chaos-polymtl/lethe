// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_linear_subequations_solver.h>
#include <solvers/vof_phase_gradient_projection.h>
#include <solvers/vof_subequations_interface.h>

template <int dim>
VOFSubequationsInterface<dim>::VOFSubequationsInterface(
  const SimulationParameters<dim> &p_simulation_parameters,
  MultiphysicsInterface<dim>      *p_multiphysics,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> &p_triangulation,
  std::shared_ptr<SimulationControl> &p_simulation_control,
  ConditionalOStream                 &p_pcout)
  : multiphysics(p_multiphysics)
  , pcout(p_pcout)
{
  if (p_simulation_parameters.multiphysics.vof_parameters.surface_tension_force
        .enable)
    {
      // Phase gradient projection
      this->active_subequations.push_back(
        VOFSubequationsID::phase_gradient_projection);
      this->subequations[VOFSubequationsID::phase_gradient_projection] =
        std::make_shared<VOFPhaseGradientProjection<
          dim,
          VOFPhaseGradientProjectionScratchData<dim>>>(this,
                                                       this->multiphysics,
                                                       p_simulation_parameters,
                                                       p_triangulation,
                                                       p_simulation_control,
                                                       this->pcout);
    }
}

template <int dim>
std::shared_ptr<PhysicsSubequationsScratchDataBase>
VOFSubequationsInterface<dim>::scratch_data_cast(
  const VOFSubequationsID  &subequation_id,
  const FiniteElement<dim> &fe_subequation,
  const Quadrature<dim>    &quadrature,
  const Mapping<dim>       &mapping,
  const FiniteElement<dim> &fe_physics)
{
  AssertThrow((std::find(this->active_subequations.begin(),
                         this->active_subequations.end(),
                         subequation_id) != this->active_subequations.end()),
              ExcInternalError());

  if (subequation_id == VOFSubequationsID::phase_gradient_projection)
    return std::make_shared<VOFPhaseGradientProjectionScratchData<dim>>(
      fe_subequation, quadrature, mapping, fe_physics);
  else // At the moment, only one option is possible. This will change with the
       // addition of other subequations to the interface.
    return std::make_shared<VOFPhaseGradientProjectionScratchData<dim>>(
      fe_subequation, quadrature, mapping, fe_physics);
}


template class VOFSubequationsInterface<2>;
template class VOFSubequationsInterface<3>;
