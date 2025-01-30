// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_algebraic_interface_reinitialization.h>
#include <solvers/vof_curvature_projection.h>
#include <solvers/vof_linear_subequations_solver.h>
#include <solvers/vof_phase_gradient_projection.h>
#include <solvers/vof_reinitialized_phase_gradient_projection.h>
#include <solvers/vof_subequations_interface.h>

template <int dim>
VOFSubequationsInterface<dim>::VOFSubequationsInterface(
  const SimulationParameters<dim> &p_simulation_parameters,
  const ConditionalOStream        &p_pcout,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> &p_triangulation,
  MultiphysicsInterface<dim> *p_multiphysics_interface)
  : multiphysics_interface(p_multiphysics_interface)
  , pcout(p_pcout)
{
  if ((p_simulation_parameters.multiphysics.vof_parameters.surface_tension_force
         .enable) ||
      (p_simulation_parameters.multiphysics.vof_parameters
         .algebraic_interface_reinitialization.enable))
    {
      // Phase gradient projection
      this->active_subequations.push_back(
        VOFSubequationsID::phase_gradient_projection);
      this->subequations[VOFSubequationsID::phase_gradient_projection] =
        std::make_shared<VOFPhaseGradientProjection<dim>>(
          p_simulation_parameters,
          this->pcout,
          p_triangulation,
          this->multiphysics_interface,
          this);
      // Phase curvature projection
      this->active_subequations.push_back(
        VOFSubequationsID::curvature_projection);
      this->subequations[VOFSubequationsID::curvature_projection] =
        std::make_shared<VOFCurvatureProjection<dim>>(
          p_simulation_parameters,
          this->pcout,
          p_triangulation,
          this->multiphysics_interface,
          this);

      if (p_simulation_parameters.multiphysics.vof_parameters
            .algebraic_interface_reinitialization.enable)
        {
          // Algebraic interface reinitialization
          this->active_subequations.push_back(
            VOFSubequationsID::algebraic_interface_reinitialization);
          this->subequations
            [VOFSubequationsID::algebraic_interface_reinitialization] =
            std::make_shared<VOFAlgebraicInterfaceReinitialization<dim>>(
              p_simulation_parameters,
              this->pcout,
              p_triangulation,
              this->multiphysics_interface,
              this);

          // Reinitialized phase gradient projection
          this->active_subequations.push_back(
            VOFSubequationsID::reinitialized_phase_gradient_projection);
          this->subequations
            [VOFSubequationsID::reinitialized_phase_gradient_projection] =
            std::make_shared<VOFReinitializedPhaseGradientProjection<dim>>(
              p_simulation_parameters,
              this->pcout,
              p_triangulation,
              this->multiphysics_interface,
              this);
        }
    }
}

template class VOFSubequationsInterface<2>;
template class VOFSubequationsInterface<3>;
