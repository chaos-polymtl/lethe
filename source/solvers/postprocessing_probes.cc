// SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#include <solvers/multiphysics_interface.h>
#include <solvers/postprocessing_probes.h>

template <int dim>
PostprocessingProbes<dim>::PostprocessingProbes(
  MultiphysicsInterface<dim> *multiphysics_interface,
  std::shared_ptr<parallel::DistributedTriangulationBase<dim>> triangulation,
  std::shared_ptr<SimulationControl> simulation_control, // for time
  const Parameters::PostProcessing<dim>::ProbingPoints
                    &probing_points_parameters,
  const std::string &output_folder)
  : multiphysics_interface(multiphysics_interface)
  , triangulation(triangulation)
  , simulation_control(simulation_control)
  , probing_points_parameters(probing_points_parameters)
  , output_folder(output_folder)
  , probe_tables(probing_points_parameters.number_of_probing_points){};

template<int dim>
void
PostprocessingProbes<dim>::postprocess_probes(Variable variable)
{

}

template class PostprocessingProbes<2>;
template class PostprocessingProbes<3>;
