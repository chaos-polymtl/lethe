// SPDX-FileCopyrightText: Copyright (c) 2024-2025 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/vof_algebraic_interface_reinitialization.h>
#include <solvers/vof_curvature_projection.h>
#include <solvers/vof_linear_subequations_solver.h>
#include <solvers/vof_phase_gradient_projection.h>
#include <solvers/vof_subequations_interface.h>

template <int dim>
void
VOFSubequationsInterface<dim>::initialize_subequations(
  const SimulationParameters<dim> &p_simulation_parameters,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  const std::shared_ptr<SimulationControl> &p_simulation_control)
{
  // Reset all subequations
  this->active_subequations.clear();
  this->subequations.clear();

  // Activate and add relevant subequations
  if ((p_simulation_parameters.multiphysics.vof_parameters.surface_tension_force
         .enable) ||
      (p_simulation_parameters.multiphysics.vof_parameters.regularization_method
         .algebraic_interface_reinitialization.enable))
    {
      // Phase gradient projection
      this->active_subequations.push_back(
        VOFSubequationsID::phase_gradient_projection);
      this->subequations[VOFSubequationsID::phase_gradient_projection] =
        std::make_shared<VOFPhaseGradientProjection<dim>>(
          p_simulation_parameters, this->pcout, p_triangulation, *this);
      // Phase curvature projection
      this->active_subequations.push_back(
        VOFSubequationsID::curvature_projection);
      this->subequations[VOFSubequationsID::curvature_projection] =
        std::make_shared<VOFCurvatureProjection<dim>>(p_simulation_parameters,
                                                      this->pcout,
                                                      p_triangulation,
                                                      *this);

      if (p_simulation_parameters.multiphysics.vof_parameters
            .regularization_method.algebraic_interface_reinitialization.enable)
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
              p_simulation_control,
              *this);
        }
    }
  reset_subequations_solutions_validity();
}


template class VOFSubequationsInterface<2>;
template class VOFSubequationsInterface<3>;
