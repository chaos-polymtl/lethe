// SPDX-FileCopyrightText: Copyright (c) 2021-2024 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/subequations_interface.h>
#include <solvers/vof_phase_gradient_projection.h>

template <int dim>
SubequationsInterface<dim>::SubequationsInterface(
  const SimulationParameters<dim>                             &sim_param,
  MultiphysicsInterface<dim>                                  *p_multiphysics,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> p_triangulation,
  std::shared_ptr<SimulationControl> p_simulation_control,
  ConditionalOStream                &p_pcout)
  : multiphysics(p_multiphysics)
  , pcout(p_pcout)
{
  if (sim_param.multiphysics.VOF)
    {
      // Phase gradient projection
      verbosity[SubequationsID::phase_gradient_projection] =
        (sim_param.non_linear_solver.at(PhysicsID::VOF).verbosity !=
           Parameters::Verbosity::quiet ||
         sim_param.linear_solver.at(PhysicsID::VOF).verbosity !=
           Parameters::Verbosity::quiet) ?
          Parameters::Verbosity::verbose :
          Parameters::Verbosity::quiet;
      active_subequations.push_back(SubequationsID::phase_gradient_projection);
      subequations[SubequationsID::phase_gradient_projection] =
        std::make_shared<VOFPhaseGradientProjection<dim>>(this->pcout,
                                                          this,
                                                          this->multiphysics,
                                                          sim_param,
                                                          p_triangulation,
                                                          p_simulation_control);
    }
}

template class SubequationsInterface<2>;
template class SubequationsInterface<3>;