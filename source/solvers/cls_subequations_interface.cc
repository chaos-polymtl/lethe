// SPDX-FileCopyrightText: Copyright (c) 2024-2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/cls_algebraic_interface_reinitialization.h>
#include <solvers/cls_curvature_projection.h>
#include <solvers/cls_linear_subequations_solver.h>
#include <solvers/cls_phase_gradient_projection.h>
#include <solvers/cls_subequations_interface.h>

template <int dim>
void
CLSSubequationsInterface<dim>::initialize_subequations(
  const SimulationParameters<dim> &p_simulation_parameters,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  const std::shared_ptr<SimulationControl> &p_simulation_control)
{
  // Reset all subequations
  this->active_subequations.clear();
  this->subequations.clear();

  // Activate and add relevant subequations
  if ((p_simulation_parameters.multiphysics.cls_parameters.surface_tension_force
         .enable) ||
      (p_simulation_parameters.multiphysics.cls_parameters.regularization_method
         .algebraic_interface_reinitialization.enable))
    {
      // Phase gradient projection
      this->active_subequations.push_back(
        CLSSubequationsID::phase_gradient_projection);
      this->subequations[CLSSubequationsID::phase_gradient_projection] =
        std::make_shared<CLSPhaseGradientProjection<dim>>(
          p_simulation_parameters, this->pcout, p_triangulation, *this);
      // Phase curvature projection
      this->active_subequations.push_back(
        CLSSubequationsID::curvature_projection);
      this->subequations[CLSSubequationsID::curvature_projection] =
        std::make_shared<CLSCurvatureProjection<dim>>(p_simulation_parameters,
                                                      this->pcout,
                                                      p_triangulation,
                                                      *this);

      if (p_simulation_parameters.multiphysics.cls_parameters
            .regularization_method.algebraic_interface_reinitialization.enable)
        {
          // Algebraic interface reinitialization
          this->active_subequations.push_back(
            CLSSubequationsID::algebraic_interface_reinitialization);
          this->subequations
            [CLSSubequationsID::algebraic_interface_reinitialization] =
            std::make_shared<CLSAlgebraicInterfaceReinitialization<dim>>(
              p_simulation_parameters,
              this->pcout,
              p_triangulation,
              p_simulation_control,
              *this);
        }
    }
  reset_subequations_solutions_validity();
}


template class CLSSubequationsInterface<2>;
template class CLSSubequationsInterface<3>;
